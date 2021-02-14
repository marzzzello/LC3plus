#!/usr/bin/env python3
#
# /******************************************************************************
# *                        ETSI TS 103 634 V1.1.1                               *
# *              Low Complexity Communication Codec Plus (LC3plus)              *
# *                                                                             *
# * Copyright licence is solely granted through ETSI Intellectual Property      *
# * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# * estoppel or otherwise.                                                      *
# ******************************************************************************/
#
# Version 1.0.4
#
# Changelog:
# 1.0.4 Better command logging
# 1.0.3 Added flag to use system sox
# 1.0.2 Reduced disk space requirement
# 1.0.1 Removed channel coder flag comparison for non-decoder tests
# 1.0.0 Initial release

import argparse
import collections
import concurrent.futures
import configparser
import datetime
import filecmp
import hashlib
import io
import itertools
import logging
import math
import os
import pathlib
import re
import shlex
import shutil
import struct
import subprocess
import sys
import uuid
import wave
import zipfile
try:
    import numpy
except ImportError:
    sys.exit('Numpy missing! Try running "pip3 install numpy".')

# convenience
Executor = concurrent.futures.ThreadPoolExecutor
TestEnv  = collections.namedtuple('TestEnv', 'profile test config test_dir, workers')
WavInfo  = collections.namedtuple('WavInfo', 'ch rate frames')

# constants
SAMPLERATES = [8000, 16000, 24000, 32000, 44100, 48000]
FRAME_SIZES = [2.5, 5.0, 10.0]
SQAM_URL    = 'https://tech.ebu.ch/docs/testmaterial/SQAM_FLAC.zip'
SQAM_SHA256 = '7d6fcd0fc42354637291792534b61bf129612f221f8efef97b62e8942a8686aa'
SOX_URL     = 'https://sourceforge.net/projects/sox/files/sox/14.4.2/sox-14.4.2-win32.zip'
SOX_SHA256  = '8072cc147cf1a3b3713b8b97d6844bb9389e211ab9e1101e432193fad6ae6662'
SOX_EXE     = pathlib.Path('SoX/sox-14.4.2/sox.exe')
INF         = float('inf')
REFERENCE_ENCODER = './LC3plus.exe -q -E -formatG192 -frame_ms {frame_size} {options} "{input}" "{output}" {bitrate}'
REFERENCE_DECODER = './LC3plus.exe -q -D -formatG192  {options} "{input}" "{output}"'

# test items
ITEM_DIR = pathlib.Path('test_items')
ITEMS = {  # start, frag, SQAM name
    'ABBA'                : ( 7,  8, '69.flac'),
    'Castanets'           : ( 0,  8, '27.flac'),
    'Eddie_Rabbitt'       : ( 0,  8, '70.flac'),
    'Female_Speech_German': ( 0,  8, '53.flac'),
    'Glockenspiel'        : ( 0, 10, '35.flac'),
    'Piano_Schubert'      : ( 0,  8, '60.flac'),
    'Violoncello'         : ( 0, 10, '10.flac'),
    'Harpsichord'         : (39,  9, '40.flac'),
    'Male_Speech_English' : ( 0,  8, '50.flac'),
}
ITEMS_PLC = ['ABBA', 'Castanets', 'Female_Speech_German', 'Harpsichord' , 'Male_Speech_English']
ITEM_LOW_PASS   = 'White_Noise_LP20'
ITEM_BAND_LIMIT = 'Female_Speech_German'

# sampling rate, band widths, bytes / frame
BAND_LIMITS = {
    48000: ([4000, 8000, 12000, 16000], 115),
    32000: ([4000, 8000, 12000], 80),
    24000: ([4000, 8000], 60),
    16000: ([4000], 40),
}
BAND_WIDTHS = {
    48000: [4000, 8000, 12000, 16000, 20000],
    32000: [4000, 8000, 12000, 16000],
    24000: [4000, 8000, 12000],
    16000: [4000, 8000],
}

# config default values
TESTS = [
    'sqam',
    'band_limiting',
    'low_pass',
    'bitrate_switching',
    'bandwidth_switching',
    'plc',
    'pc',
    'ep_correctable',
    'ep_non_correctable',
    'ep_mode_switching',
    'ep_combined',
    'ep_combined_nc',
]
TEST_MODES = ['encode', 'encdec', 'decode']
DEFAULTS_GLOBAL = {
    'option_bandwidth' : '',
    'option_ep_debug'  : '',
    'option_ep_mode'   : '',
    'option_plc_mode'  : '',
    'peaq_bin'         : '',
    'peaq_odg_regex'   : '',
    'reference_decoder': REFERENCE_DECODER,
    'reference_encoder': REFERENCE_ENCODER,
}
DEFAULTS_TEST = {'configs': []}
for test in TESTS:
    DEFAULTS_TEST['test_' + test] = False
for test, mode in itertools.product(TESTS, TEST_MODES):
    DEFAULTS_TEST['{}_{}_eng_threshold'.format(test, mode)] = 70
    DEFAULTS_TEST['{}_{}_mld_threshold'.format(test, mode)] = 4
    DEFAULTS_TEST['{}_{}_odg_threshold'.format(test, mode)] = 0.06
    DEFAULTS_TEST['{}_{}_rms_threshold'.format(test, mode)] = 14
    DEFAULTS_TEST['{}_{}_metric'.format(test, mode)]        = 'rms'
METRIC_DEFAULTS = {
    'low_pass_encode_metric'          : 'eng',
    'low_pass_encdec_metric'          : 'eng',
    'plc_decode_metric'               : 'mld',
    'pc_decode_metric'                : 'mld',
    'ep_non_correctable_decode_metric': 'mld',
    'ep_non_correctable_encdec_metric': 'mld',
    'ep_combined_nc_decode_metric'    : 'mld',
    'ep_combined_nc_encdec_metric'    : 'mld',
}
DEFAULTS_TEST.update(METRIC_DEFAULTS)

# html output stuff
LABEL = {
    'sqam'               : 'SQAM',
    'band_limiting'      : 'Band Limitation',
    'low_pass'           : 'Low Pass',
    'bitrate_switching'  : 'Bitrate Switching',
    'bandwidth_switching': 'Bandwidth Switching',
    'plc'                : 'Packet Loss Concealment',
    'pc'                 : 'Partial Concealment',
    'ep_correctable'     : 'Channel Coder for Correctable Frames',
    'ep_non_correctable' : 'Channel Coder for Non-Correctable Frames',
    'ep_mode_switching'  : 'Error protection mode switching',
    'ep_combined'        : 'Combined Channel Coder for Correctable Frames',
    'ep_combined_nc'     : 'Combined Channel Coder for Non-Correctable Frames',
}
HEADER_ALL = ['Mode', 'Item', 'Frame Size', 'Samplerate', 'Bitrate']
HEADER_EP  = ['EP Mode']
HEADER_EPD = ['BFI', 'EPMR', 'ER']
HEADER_BL  = ['Bandwidth']
HEADER_METRIC = {
    'rms': ['Max. Abs. Diff', 'RMS [dB]', 'RMS Reached [bits]'],
    'odg': ['ODG<sub>ref</sub>', 'Δ<sub>ODG</sub>'],
    'mld': ['MLD'],
    'eng': ['E<sub>diff</sub> [dB]'],
    None : []
}
HTML_HEAD = ('<!DOCTYPE html><head><meta charset="UTF-8"><title>{title} Report</title>'
             '<style>{style}</style></head>\n<body><h2>Conformance for "{title}" {state}!</h2>\n')
HTML_DIV = '<div><h3>{label} - {percent}%</h3>\n'
STYLE = ('body {font-family:sans-serif; color:#f8f8f2; background-color:#272822; font-size:80%} div {border:1px solid '
         '#8f908a; border-radius:4px; overflow:hidden; display:table; margin-left:30px; margin-bottom:30px} h2 {text-a'
         'lign:left; margin-left:30px} h3 {text-align:left; margin:4px} table {border-spacing:0px; width:100%} th {pad'
         'ding:4px; border-top:1px solid #8f908a} td {padding:4px} tr:nth-child(even) {background-color:rgba(255,255,2'
         '55,0.1)} td.pass {background-color:rgba(0,192,255,0.4)} td.fail {background-color:rgba(255,0,0,0.4)} td.warn'
         '{background-color:rgba(214,137,16,0.4)}')


# convenience wrapper for os.makedirs
def makedirs(path):
    os.makedirs(str(path), exist_ok=True)
    return path


# convenience wrapper for shutil.rmtree
def removedir(path):
    shutil.rmtree(str(path), ignore_errors=True)


# Run command and return output. cmd can be string or list. Commands with .exe suffix are automatically
# called with wine unless wine=False. Set unicode=False to get binary output. Set hard_fail=False to
# to ignore nonzero return codes.
def call(cmd, wine=True, unicode=True, hard_fail=True, log_output=True):
    if isinstance(cmd, str):
        cmd = [x for x in shlex.split(cmd) if x]
    if sys.platform != 'cygwin' and wine and cmd[0].lower().endswith('.exe'):
        cmd = ['wine'] + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=unicode)
    out = p.communicate()[0] or (b'', '')[unicode]
    quoted_cmd = ' '.join(map(shlex.quote, cmd))
    logging.debug(quoted_cmd)
    if unicode and log_output:
        logging.debug(out)
    if hard_fail and p.returncode != 0:
        raise OSError(quoted_cmd + ' failed!')
    return out


# return url as as bytes object, validate against hash
def download(url, sha256=None):
    try:
        buf = call('curl --silent -L "{}"'.format(url), unicode=False)
    except OSError:
        sys.exit('Failed to download {}!'.format(url))
    if sha256 and hashlib.sha256(buf).hexdigest() != sha256:
        sys.exit('Failed to validate hash for {}!'.format(url))
    return buf


def download_sox():
    if not is_file(SOX_EXE):
        print('Downloading SoX ...')
        buf = download(SOX_URL, SOX_SHA256)
        zipfile.ZipFile(io.BytesIO(buf)).extractall(str(SOX_EXE.parent.parent))
        if sys.platform == 'cygwin':
            call('chmod -R +x "{}"'.format(SOX_EXE.parent))


def exe_exists(exe, wine=False):
    try:
        out = call(exe, wine=wine, hard_fail=False)
    except OSError:
        return False
    return not (wine and out.startswith('wine: ')) # detect wine: cannot find


def check_system(args, globvars):
    if sys.platform == 'win32':
        sys.exit('This script must run under cygwin')
    tools = ['curl', 'gcc', 'make']
    if sys.platform != 'cygwin':
        tools += ['wine']
    if args.system_sox:
        tools += ['sox']
    for tool in tools:
        if not exe_exists(tool):
            sys.exit('Failed to find {} executable'.format(tool))
    if globvars['peaq_bin'] and not exe_exists(globvars['peaq_bin'], wine=True):
        sys.exit('Failed to find PEAQ executable. Please adjust config file')
    if not exe_exists(globvars['encoder'], wine=True):
        sys.exit('Failed to find LC3 encoder executable. Please adjust config file.')
    if not exe_exists(globvars['decoder'], wine=True):
        sys.exit('Failed to find LC3 decoder executable. Please adjust config file.')


# search s for expr and return first match or exit
def regex_search(expr, s):
    if not re.search(expr, s):
        sys.exit('No match for regular expression "{}"!'.format(expr))
    return re.search(expr, s).group(1)


# calculates the max xcorr of the two vectors within the length of the longer one
def align_vec(x1, x2):
    # trims longer vector at position of max xcorr and returns both
    assert len(x1) >= len(x2)
    res = []
    d = len(x1) - len(x2)
    # normalize to max of int16
    a = numpy.float32(x1) / 32767
    b = numpy.float32(x2) / 32767
    # a is longer than b
    for i in range(d + 1):
        xx = numpy.dot(a[i:len(b) + i], b)
        res.append(xx)
    lag = numpy.array(res).argmax()
    # trim longer vector
    x1 = x1[lag:lag + len(x2)]
    return x1, x2


# convert byte objects to signed int16
def byte_to_float(b, frames, channels):
    return struct.unpack("%ih" % (frames * channels), b)


# trim longer file to length of shorter file at max xcorr position and overwrite longer one
def align_files(file_1, file_2):
    logging.debug('align_files: %s %s', file_1, file_2)
    file_1, file_2 = str(file_1), str(file_2)
    # read in audio files and to max of int16
    if wav_info(file_1).frames < wav_info(file_2).frames:
        file_1, file_2 = file_2, file_1
    with wave.open(file_1, 'rb') as wf1, wave.open(file_2, 'rb') as wf2:
        b1 = wf1.readframes(wf1.getnframes())
        b2 = wf2.readframes(wf2.getnframes())
        x1 = byte_to_float(b1, wf1.getnframes(), wf1.getnchannels())
        x2 = byte_to_float(b2, wf2.getnframes(), wf2.getnchannels())
        par1 = wf1.getparams()
        par2 = wf2.getparams()
    # measure cross correlation -> delay between files and return trimmed vectors
    y1, y2 = align_vec(x1, x2)
    # overwrite files
    with wave.open(file_1, 'wb') as wf1, wave.open(file_2, 'wb') as wf2:
        wf1.setparams(par1)
        wf2.setparams(par2)
        wf1.setnframes(len(y1))
        wf2.setnframes(len(y2))
        b1 = struct.pack("%ih" % len(y1), *y1)
        b2 = struct.pack("%ih" % len(y2), *y2)
        wf1.writeframes(b1)
        wf2.writeframes(b2)


def build_tools():
    call('make -C tools')


# return info about wav file
def wav_info(path):
    wav = wave.open(str(path))
    return WavInfo(wav.getnchannels(), wav.getframerate(), wav.getnframes())


# call sox with args in repeatable mode, lazy skips execution if output already exists
def sox(*args, lazy=False):
    wavs = [x for x in map(str, args) if x.endswith('.wav')]
    if not (lazy and is_file(wavs[-1])): # last .wav is assumed to be output
        call('{} -R {}'.format(SOX_EXE, ' '.join(map(str, args))))


def trim(input, output, start, end, lazy=False):
    if not (lazy and is_file(output)):
        tmp = output.parent / uuid_file('trim_', '.wav')
        sox(input, tmp, 'trim', start, end)
        wi = wav_info(tmp)
        sox(tmp, output, 'fade', 0.5, wi.frames / wi.rate, 0.7)
        tmp.unlink()


# resample wav using sox
def resample(input, output, rate, lazy=False):
    sox(input, output, 'rate -vs', rate, lazy=lazy)


# apply lowpass filter using sox
def low_pass(input, output, fc, lazy=False):
    sox(input, output, 'sinc -{}'.format(fc), lazy=lazy)


# generate switching file with unique name, returns path
def generate_switching_file(env, *values):
    path   = env.test_dir / uuid_file('swf_', '.dat')
    layers = ','.join(map(str, sorted(values)))
    cmd    = 'tools/gen-rate-profile.exe -layers {} {}'
    call(cmd.format(layers, path), log_output=False)
    return path


# compares binary equality of files
def compare_bin(file1, file2):
    logging.debug('compare_bin: %s %s', file1, file2)
    return filecmp.cmp(str(file1), str(file2))


# copy file from src to dst
def copy_file(src, dst):
    logging.debug('copy_file: %s %s', src, dst)
    shutil.copy(str(src), str(dst))


# generate unique file name with extension
def uuid_file(prefix='', suffix=''):
    return prefix + str(uuid.uuid4()) + suffix


# like str() but with special case for float
def fstr(x):
    return '{:.3g}'.format(x) if type(x) == float else str(x)


# like str() but with special case for list
def lstr(x):
    return '-'.join(map(str, x)) if type(x) in (list, tuple) else str(x)


# returns true if path is a file
def is_file(path):
    return os.path.isfile(str(path))


# calculate bitrate from bytes per frame
def get_bitrate(bytes_per_frame, frame_size):
    return int(bytes_per_frame * 8000 / frame_size)


# apply func to list of argumets,
def thread_executor(func, args, workers):
    list(ThreadPoolExecutor(workers).map(lambda x: func(*x), args)) # list() to collect futures


def prepare_items(workers):
    sqam_dir = pathlib.Path('SQAM')
    item_dir = makedirs(ITEM_DIR)
    if not sqam_dir.exists():
        print('Downloading test items ...')
        buf = download(SQAM_URL, SQAM_SHA256)
        zipfile.ZipFile(io.BytesIO(buf)).extractall(str(sqam_dir))

    print('Preparing test items ...')
    # trim items
    with Executor(workers) as ex:
        for item, (st, fr, flac) in ITEMS.items():
            infile  = sqam_dir / flac
            outfile = item_dir / (item + '.wav')
            ex.submit(trim, infile, outfile, st, fr, lazy=True)
    # resample items
    with Executor(workers) as ex:
        for item, sr in itertools.product(ITEMS, SAMPLERATES):
            infile  = item_dir / (item + '.wav')
            outfile = item_dir / '{}_{}_2ch.wav'.format(item, sr)
            ex.submit(resample, infile, outfile, sr, lazy=True)
    with Executor(workers) as ex:
        # 20 kHz lowpass
        for item, sr in itertools.product(ITEMS, SAMPLERATES):
            if sr >= 44100:
                infile  = item_dir / '{}_{}_2ch.wav'.format(item, sr)
                outfile = item_dir / '{}_{}_2ch_lp20.wav'.format(item, sr)
                ex.submit(low_pass, infile, outfile, 20000, lazy=True)
        # band limit
        for sr, (bws, _) in BAND_LIMITS.items():
            for bw in bws:
                infile  = item_dir / '{}_{}_2ch.wav'.format(ITEM_BAND_LIMIT, sr)
                outfile = item_dir / '{}_{}_2ch_bw{}.wav'.format(ITEM_BAND_LIMIT, sr, bw)
                ex.submit(low_pass, infile, outfile, bw, lazy=True)
        # LP20 item with 4 seconds of white noise above 20kHz
        outfile = item_dir / (ITEM_LOW_PASS + '_48000_1ch.wav')
        synth   = 'synth 4 white fir etc/hp_fir_coef.txt'
        ex.submit(sox, '-n -r 48000 -c 1 -b 16', outfile, synth, lazy=True)
    # create 1ch items
    with Executor(workers) as ex:
        for path in os.listdir(str(item_dir)):
            if '_2ch' in path:
                infile  = item_dir / path
                outfile = item_dir / path.replace('_2ch', '_1ch')
                ex.submit(sox, infile, outfile, 'remix -', lazy=True)


def parse_config(path):
    def strip_comment(line):
        return line.split('#', 1)[0].strip()

    def split_list(line):
        return [x.strip() for x in strip_comment(line).split(',')]

    def parse_conf_line(line):
        mode, fs, sr, br = split_list(line)
        fs, sr = float(fs), int(sr)
        if ':' in br:
            br_start, br_step, br_stop = map(int, br.split(':'))
            br = list(range(br_start, br_stop + 1, br_step))
        else:
            br = [int(br)]
        if fs not in FRAME_SIZES:
            sys.exit('Unsupported frame size: {}!'.format(line))
        if sr not in SAMPLERATES:
            sys.exit('Unsupported sampling rate: {}!'.format(line))
        if min(br) < 16000 or max(br) > 320000:
            sys.exit('Invalid bitrate: {}!'.format(line))
        return mode, fs, sr, br

    def parse_bool(val):
        if val not in ('0', '1'):
            raise ValueError
        return val == '1'

    if not is_file(path):
        sys.exit('No such file: ' + path)

    globals_required = ['enabled_tests', 'encoder', 'decoder']
    globals_all      = list(DEFAULTS_GLOBAL) + globals_required
    test_keys = ['test_' + t for t in TESTS]
    globels   = DEFAULTS_GLOBAL.copy()
    configs   = {}

    try:
        parser = configparser.ConfigParser()
        parser.read(path)
        # parse global section
        for key in parser['globals']:
            globels[key] = strip_comment(parser['globals'][key])
            if key not in globals_all:
                sys.exit('Unknown key "{}" in config'.format(key))
        globels['enabled_tests'] = split_list(parser['globals']['enabled_tests'])
        # trigger KeyError for required keys
        map(lambda key: globels[key], globals_required)
        # parse test sections
        for profile in globels['enabled_tests']:
            configs[profile] = {**globels, **DEFAULTS_TEST}
            for key in parser[profile]:
                try:
                    val = strip_comment(parser[profile][key])
                    if key in test_keys:
                        configs[profile][key] = parse_bool(val)
                    elif key == 'configs':
                        lines = parser[profile][key].splitlines()
                        configs[profile][key] = [parse_conf_line(l) for l in lines]
                    elif key.endswith('_threshold') and key in DEFAULTS_TEST:
                        configs[profile][key] = float(val)
                    elif key.endswith('_metric') and key in DEFAULTS_TEST:
                        if val not in ('rms', 'odg', 'mld', 'eng'):
                            raise ValueError
                        configs[profile][key] = val
                    else:
                        sys.exit('Unknown key "{}" in config'.format(key))
                except ValueError:
                    sys.exit('Invalid value in config: {} = {}'.format(key, parser[profile][key]))
    except KeyError as e:
        sys.exit('Missing "{}" in config'.format(e.args[0]))
    except configparser.DuplicateOptionError as e:
        sys.exit('Duplicate key "{}" in config'.format(e.args[1]))

    return globels, configs


# splits up files into channels, yields channel files
def split_channels(env, *files):
    channels = wav_info(files[0]).ch
    if channels == 1:
        yield files
    else:
        for ch in range(1, channels + 1):
            tmp_files = []
            for f in files:
                tmp = env.test_dir / uuid_file('split_', '.wav')
                sox(f, tmp, 'remix', ch)
                tmp_files.append(tmp)
            yield tmp_files


def run_rms(env, file1, file2, threshold):
    rms, diff, bits = -INF, 0, 24
    for split1, split2 in split_channels(env, file1, file2):
        tmp1 = env.test_dir / uuid_file('rms_', '.wav')
        tmp2 = env.test_dir / uuid_file('rms_', '.wav')
        copy_file(split1, tmp1)
        copy_file(split2, tmp2)
        align_files(tmp1, tmp2)
        out = call('tools/rms {} {} {}'.format(tmp1, tmp2, threshold))
        diff_samp = int(regex_search(r'different samples\s+: (\d+)', out))
        if diff_samp != 0:
            rms  = max(rms, float(regex_search(r'Overall RMS value\s+: (\S+) dB ---', out)))
            diff = max(diff, float(regex_search(r'Maximum difference\s+: (\S+) ---', out)))
            bits = min(bits, int(regex_search(r'RMS criteria\s+: (\d+) bit', out)))
    return rms, diff, bits


def run_peaq(env, reference, test):
    odg = 5
    for split1, split2 in split_channels(env, reference, test):
        ref = env.test_dir / uuid_file('odg_', '.wav')
        tst = env.test_dir / uuid_file('odg_', '.wav')
        resample(split1, ref, 48000)
        resample(split2, tst, 48000)
        align_files(ref, tst)
        out = call(env.config['peaq_bin'].format(reference=ref, test=tst))
        odg = min(odg, float(regex_search(env.config['peaq_odg_regex'], out)))
    return odg


def run_mld(env, reference, test):
    mld = 0
    for split1, split2 in split_channels(env, reference, test):
        ref = env.test_dir / uuid_file('mld_', '.wav')
        tst = env.test_dir / uuid_file('mld_', '.wav')
        resample(split1, ref, 48000)
        resample(split2, tst, 48000)
        align_files(ref, tst)
        out = call('tools/mld -d {} {}'.format(ref, tst))
        mld = max(mld, float(regex_search(r'maximum loudness difference:\s*(\S+)', out)))
    return mld


# calculate energy difference of two wavs
def energy_diff(env, file1, file2):
    logging.debug('energy_diff: %s %s', file1, file2)
    tmp1 = str(env.test_dir / uuid_file('eng_', '.wav'))
    tmp2 = str(env.test_dir / uuid_file('eng_', '.wav'))
    copy_file(file1, tmp1)
    copy_file(file2, tmp2)
    align_files(tmp1, tmp2)
    with wave.open(tmp1, 'rb') as ref_wf, wave.open(tmp2, 'rb') as tst_wf:
        bytes_ref = ref_wf.readframes(ref_wf.getnframes())
        bytes_tst = tst_wf.readframes(tst_wf.getnframes())
        ref = byte_to_float(bytes_ref, ref_wf.getnframes(), ref_wf.getnchannels())
        tst = byte_to_float(bytes_tst, tst_wf.getnframes(), tst_wf.getnchannels())
        eng_diff = sum(numpy.square(numpy.subtract(ref, tst)))
        eng_diff = math.log10(eng_diff) if eng_diff != 0 else -INF
        return eng_diff


# compare output wavs by metric rms, odg, mld, eng
def compare_wav(env, mode, infile, file_ref, file_tst):
    mkey   = '{}_{}_metric'.format(env.test, mode)
    metric = env.config[mkey]
    tkey   = '{}_{}_{}_threshold'.format(env.test, mode, metric)
    thresh = env.config[tkey]

    if metric == 'rms':
        rms, diff, bits = run_rms(env, file_ref, file_tst, thresh)
        rms_thr  = 20 * math.log10(2 ** (-thresh + 1) / 12 ** 0.5)
        diff_thr = 1 / 2 ** (thresh - 3)
        ok_rms   = rms <= rms_thr
        ok_diff  = diff <= diff_thr
        ok_bits  = bits >= thresh
        ok       = ok_rms and ok_diff
        values   = [(diff,  ('fail', 'pass')[ok_diff], diff_thr),
                    (rms, ('fail', 'pass')[ok_rms], rms_thr),
                    (bits, ('warn', 'none')[ok_bits], thresh)]
    if metric == 'odg':
        odg_ref  = run_peaq(env, infile, file_ref)
        odg_tst  = run_peaq(env, infile, file_tst)
        odg_diff = abs(odg_ref - odg_tst)
        ok       = odg_diff <= thresh
        values   = [(odg_ref, '', None),
                    (odg_diff, ('fail', 'pass')[ok], thresh)]
    if metric == 'mld':
        mld    = run_mld(env, file_ref, file_tst)
        ok     = mld <= thresh
        values = [(mld, ('fail', 'pass')[ok], thresh)]
    if metric == 'eng':
        d_eng  = energy_diff(env, file_ref, file_tst)
        ok     = d_eng <= thresh
        values = [(d_eng, ('fail', 'pass')[ok], thresh)]

    return ok, values


# compare outout files of ep debug flag
def compare_errors(env, file_ref, file_tst):
    ok_all, values = True, []
    for ext in ['.bfi', '.epmr', '.error_report']:
        ok = compare_bin(str(file_ref) + ext, str(file_tst) + ext)
        values += [(('bad', 'ok')[ok], ('fail', 'pass')[ok], None)]
    return ok_all, values


# ensure inputs exist and nothig is overwritten we're not overwriting anythong
def check_io_files(input, output):
    if not is_file(input):
        raise FileNotFoundError(input)
    if is_file(output):
        raise FileExistsError(output)


def encode_reference(env, input, output, frame_size, bitrate, bandwidth=None, ep_mode=0):
    check_io_files(input, output)
    cmd = env.config['reference_encoder']
    opt = []
    if bandwidth:
        opt += ['-bandwidth', bandwidth]
    if ep_mode:
        opt += ['-epmode', ep_mode]
    options = ' '.join(map(str, opt))
    call(cmd.format(input=input, output=output, frame_size=frame_size, bitrate=bitrate, options=options))


def encode_test(env, input, output, frame_size, bitrate, bandwidth=None, ep_mode=0):
    check_io_files(input, output)
    cmd = env.config['encoder']
    opt = []
    if bandwidth:
        opt += [env.config['option_bandwidth'].format(arg=bandwidth)]
    if ep_mode:
        opt += [env.config['option_ep_mode'].format(arg=ep_mode)]
    options = ' '.join(opt)
    call(cmd.format(input=input, output=output, frame_size=frame_size, bitrate=bitrate, options=options))


def decode_reference(env, input, output, error_file=None):
    check_io_files(input, output)
    cmd = env.config['reference_decoder']
    opt = []
    if error_file:
        opt += ['-ep_dbg', error_file]
    options = ' '.join(map(str, opt))
    call(cmd.format(input=input, output=output, options=options))


def decode_test(env, input, output, error_file=None):
    check_io_files(input, output)
    cmd = env.config['decoder']
    opt = []
    if error_file:
        opt += [env.config['option_ep_debug'].format(arg=error_file)]
    options = ' '.join(map(str, opt))
    call(cmd.format(input=input, output=output, options=options))


def apply_error_pattern(env, input, output, mode, pattern):
    assert mode in ('fer', 'ber', 'flip')
    check_io_files(input, output)
    if mode == 'fer':
        cmd = 'tools/eid-xor.exe -vbr -bs g192 -ep byte -fer {} {} {}'
        call(cmd.format(input, pattern, output))
    if mode == 'ber':
        cmd = 'tools/eid-xor.exe -vbr -bs g192 -ep byte -ber {} {} {}'
        call(cmd.format(input, pattern, output))
    if mode == 'flip':
        cmd = 'tools/flipG192 {} {} {} {} 1911 0'
        flips, frames = pattern
        call(cmd.format(input, output, flips, frames))
    # copy the config file of g192 bitstreams
    if is_file(str(input) + '.cfg'):
        copy_file(str(input) + '.cfg', str(output) + '.cfg')


# create file names for test
def make_files(env, files, *args):
    protoyp = '_'.join(map(lstr, args)) + '_'
    return tuple(env.test_dir / (protoyp + f) for f in files)


# permutate test configs
def sqam_configs(config, items, ch=1, lp20=True, modes=None):
    for mode, fs, sr, brs in config['configs']:
        if modes and mode not in modes:
            continue
        for item, br in itertools.product(items, brs):
            suffix = '_lp20' if lp20 and sr in (44100, 48000) else ''
            infile = ITEM_DIR / '{}_{}_{}ch{}.wav'.format(item, sr, ch, suffix)
            yield mode, item, fs, sr, br, infile


# apply test func to list of tests, multithreadded
def test_executor(env, func, tests):
    ex = Executor(env.workers)
    return list(ex.map(lambda args: func(*args), tests))


# process a single test item
# performs encoding, erroring, decoding and evaluation
# returns tuble of bool, list (pass condition, metric values)
# mode: encode/decode/encdec, fs: frame size, sr: sampling rate, br: bitrate
def process_item(env, mode, item, fs, sr, br, infile, bandwidth=None, ep_mode=0, error_mode=None, error_pattern=None):
    bw_name = 'swf' if is_file(bandwidth) else bandwidth
    ep_name = '1-4' if is_file(ep_mode) else ep_mode
    fmt     = '  {} {:20} {:3g} ms {:5} Hz {:>6} bit/s ep:{}'
    print(fmt.format(mode, item, fs, sr, lstr(br), ep_name))

    file_names = ['r.g192', 't.g192', 're.g192', 'te.g192', 'r.wav', 't.wav', 'rd', 'td']
    file_tuple = make_files(env, file_names, mode, infile.stem, fs, sr, br, bw_name, ep_name)
    ref_bin, tst_bin, ref_err, tst_err, ref_wav, tst_wav, ref_dbg, tst_dbg = file_tuple
    # evaluate channel coder output only for decode_ep_* tests
    if not (env.test.startswith('ep_') and mode == 'decode'):
        err_ok, err_val, ref_dbg, tst_dbg = True, [], None, None

    try:
        # generate bitratre switching file if needed
        if type(br) in (list, tuple):
            br = generate_switching_file(env, *br)
        # encode
        encode_reference(env, infile, ref_bin, fs, br, bandwidth=bandwidth, ep_mode=ep_mode)
        if mode in ('encode', 'encdec'):
            encode_test(env, infile, tst_bin, fs, br, bandwidth=bandwidth, ep_mode=ep_mode)
        # apply errors
        if error_mode:
            apply_error_pattern(env, ref_bin, ref_err, error_mode, error_pattern)
            ref_bin = ref_err
            if mode in ('encode', 'encdec'):
                apply_error_pattern(env, tst_bin, tst_err, error_mode, error_pattern)
                tst_bin = tst_err
        # decode
        decode_reference(env, ref_bin, ref_wav, error_file=ref_dbg)
        if mode == 'encode':
            decode_reference(env, tst_bin, tst_wav, error_file=tst_dbg)
        if mode == 'encdec':
            decode_test(env, tst_bin, tst_wav, error_file=tst_dbg)
        if mode == 'decode':
            decode_test(env, ref_bin, tst_wav, error_file=tst_dbg)
        # compare outputs
        ok, val = compare_wav(env, mode, infile, ref_wav, tst_wav)
        if ref_dbg and tst_dbg:
            err_ok, err_val = compare_errors(env, ref_dbg, tst_dbg)
        return ok and err_ok, err_val + val

    except (OSError, FileNotFoundError, FileExistsError, KeyError) as e:
        logging.error('process_item: %s: %s', type(e).__name__, str(e))
        return False, []


def test_sqam(env):
    print('Testing SQAM ...')
    def func(mode, item, fs, sr, br, infile):
        ok, values = process_item(env, mode, item, fs, sr, br, infile)
        return [ok, mode, item, fs, sr, br] + values

    tests = sqam_configs(env.config, ITEMS)
    return test_executor(env, func, tests)


def test_band_limiting(env):
    print('Testing band limitation ...')
    def func(mode, item, fs, sr, br, bw):
        infile = ITEM_DIR / '{}_{}_1ch_bw{}.wav'.format(item, sr, bw)
        ok, values = process_item(env, mode, item, fs, sr, br, infile)
        return [ok, mode, item, fs, sr, br, bw] + values

    tests = set()
    for mode, fs, sr, _ in env.config['configs']:
        if sr >= 16000:
            bw_limits, frame_bytes = BAND_LIMITS[sr]
            br = get_bitrate(frame_bytes, fs)
            for bw in bw_limits:
                tests.add((mode, ITEM_BAND_LIMIT, fs, sr, br, bw))
    return test_executor(env, func, list(tests))


def test_low_pass(env):
    print('Testing low pass ...')
    def func(mode, item, fs, sr, br, infile):
        ok, values = process_item(env, mode, item, fs, sr, br, infile)
        return [ok, mode, item, fs, sr, br] + values

    items = [ITEM_LOW_PASS]
    modes = ['encode', 'encdec']
    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, items , modes=modes, lp20=False):
        if sr >= 44100:
            tests.append((mode, item, fs, sr, br, infile))
    return test_executor(env, func, tests)


def test_bitrate_switching(env):
    print('Testing bitrate switching ...')
    def func(mode, item, fs, sr, bitrates):
        br     = (int(160000 / fs), max(bitrates))
        infile = ITEM_DIR / '{}_{}_1ch.wav'.format(item, sr)
        ok, values = process_item(env, mode, item, fs, sr, br, infile)
        return [ok, mode, item, fs, sr, br] + values

    tests = []
    for mode, fs, sr, bitrates in env.config['configs']:
        for item in ITEMS:
            tests.append((mode, item, fs, sr, bitrates))
    return test_executor(env, func, tests)


def test_bandwidth_switching(env):
    print('Testing bandwidth switching ...')
    def func(mode, item, fs, sr, br, infile):
        bwf = generate_switching_file(env, *BAND_WIDTHS[sr])
        ok, values = process_item(env, mode, item, fs, sr, br, infile, bwf)
        return [ok, mode, item, fs, sr, br] + values

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS):
        if sr >= 16000:
            tests.append((mode, item, fs, sr, br, infile))
    return test_executor(env, func, tests)


def test_plc(env):
    print('Testing packet loss concealment ...')
    def func(mode, item, fs, sr, br, infile):
        pattern = 'etc/plc_fer_eid.dat'
        ok, values = process_item(env, mode, item, fs, sr, br, infile, None, 0, 'fer', pattern)
        return [ok, mode, item, fs, sr, br] + values

    tests = sqam_configs(env.config, ITEMS_PLC, modes=['decode'])
    return test_executor(env, func, tests)


def test_pc(env):
    print('Testing partial concealment ...')
    def func(mode, item, fs, sr, br, infile):
        pattern = 'etc/pc_ber_3percent.dat'
        ok, values = process_item(env, mode, item, fs, sr, br, infile, None, 4, 'ber', pattern)
        return [ok, mode, item, fs, sr, br] + values

    tests = sqam_configs(env.config, ITEMS_PLC, modes=['decode'])
    return test_executor(env, func, tests)


def test_ep_correctable(env):
    print('Testing channel coder for correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = (ep_mode - 1, 50)
        ok, values = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'flip', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS):
        for ep_mode in [1, 2, 3, 4]:
            tests.append((mode, item, fs, sr, br, infile, ep_mode))
    return test_executor(env, func, tests)


def test_ep_non_correctable(env):
    print('Testing channel coder for non-correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = (int(br * ep_mode * fs / 24000), 50)
        ok, values = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'flip', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS):
        for ep_mode in [1, 2, 3, 4]:
            tests.append((mode, item, fs, sr, br, infile, ep_mode))
    return test_executor(env, func, tests)


def test_ep_mode_switching(env):
    print('Testing ep-mode switching ...')
    ep_mode = generate_switching_file(env, 100, 200, 300, 400)
    def func(mode, item, fs, sr, br, infile):
        ok, values = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, None, None)
        return [ok, mode, item, fs, sr, br, '1-4'] + values

    tests = sqam_configs(env.config, ITEMS)
    return test_executor(env, func, tests)


def test_ep_combined(env):
    print('Testing combined channel coder for correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = (ep_mode - 1, 50)
        ok, values = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'flip', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS, ch=2):
        if br <= 128000:
            for ep_mode in [1, 2, 3, 4]:
                tests.append((mode, item, fs, sr, br, infile, ep_mode))

    return test_executor(env, func, tests)


def test_ep_combined_nc(env):
    print('Testing combined channel coder for non-correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = (int(br * ep_mode * fs / 24000), 50)
        ok, values = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'flip', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS, ch=2):
        if br <= 128000:
            for ep_mode in [1, 2, 3, 4]:
                tests.append((mode, item, fs, sr, br, infile, ep_mode))
    return test_executor(env, func, tests)


def pass_ratio(results):
    num_passed = sum(ok for ok, *_ in results)
    return num_passed / len(results) if results else 1.0


def profile_passed(all_results):
    flat_results = itertools.chain(*all_results.values())
    return all(ok for ok, *_ in flat_results)


def gen_td(value):
    if type(value) in (tuple, list) and len(value) == 3:
        value, clazz, thresh = value
        clazz  = ' class={}'.format(clazz) if clazz else ''
        thresh = ' ({})'.format(fstr(thresh)) if thresh != None else ''
        return '<td{}>{}{}</td>'.format(clazz, fstr(value), thresh)
    else:
        return '<td>{}</td>'.format(lstr(value))


def gen_table(test, mode, config, results):
    mkey   = '{}_{}_metric'.format(test, mode)
    metric = config[mkey]
    header = HEADER_ALL.copy()
    if test == 'band_limiting':
        header += HEADER_BL
    if test.startswith('ep_'):
        header += HEADER_EP
        if mode == 'decode':
            header += HEADER_EPD
    header += HEADER_METRIC[metric]
    buf = '<table>\n<tr>'
    buf += ''.join('<th>{}</th>'.format(x) for x in header)
    buf += '</tr>\n'
    for values in results:
        buf += '<tr>'
        buf += ''.join(map(gen_td, values[1:]))
        buf += '</tr>\n'
    return buf + '</table>\n'


def gen_div(test, config, results):
    percent = round(100 * pass_ratio(results))
    buf     = HTML_DIV.format(label=LABEL[test], percent=percent)
    for mode in TEST_MODES:
        mode_results = [r for r in results if r[1] == mode]
        if mode_results:
            buf += gen_table(test, mode, config, mode_results)
    return buf + '</div>\n'


def save_html(profile, config, all_results, html_file):
    ok    = profile_passed(all_results)
    state = ('failed', 'passed')[ok]
    buf   = HTML_HEAD.format(title=profile, style=STYLE, state=state)
    for test in TESTS:
        if test in all_results:
            buf += gen_div(test, config, all_results[test])
    buf += '</body>\n'

    with open(html_file, 'w') as f:
        f.write(buf)


def main(args):
    args.workers = min(max(args.workers, 1), os.cpu_count())
    time_stamp   = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    log_file     = 'lc3_conformance_{}.log'.format(time_stamp)
    log_handlers = [logging.FileHandler(log_file)]
    if args.verbose:
        log_handlers += [logging.StreamHandler(sys.stdout)]
    logging.basicConfig(level=logging.DEBUG, handlers=log_handlers)
    work_dir = makedirs(pathlib.Path('lc3_conformance_' + time_stamp))

    try:
        all_passed = True
        globels, configs = parse_config(args.config)
        profiles = globels['enabled_tests']
        check_system(args, globels)
        build_tools()
        if not args.system_sox:
            download_sox()
        prepare_items(args.workers)

        for profile in profiles:
            print('Running tests for "{}" ...'.format(profile))
            config      = configs[profile]
            all_results = {}
            for test in TESTS:
                test_test = 'test_' + test
                if config[test_test]:
                    test_dir    = makedirs(work_dir / profile / test)
                    test_env    = TestEnv(profile, test, config, test_dir, args.workers)
                    test_func   = globals()[test_test]
                    test_result = test_func(test_env)
                    if not test_result:
                        print('{} in "{}" is enabled with no suitable configuration!'.format(test, profile))
                    all_results[test] = test_result
                    if not args.keep_files:
                        removedir(work_dir / profile / test)

            if all_results:
                all_passed = all_passed and profile_passed(all_results)
                html_file  = '{}_{}.html'.format(profile, time_stamp)
                print('Saving results ...')
                save_html(profile, config, all_results, html_file)
            else:
                print('No tests in "{}" were enabled!'.format(profile))

        print('\nLogfile:', log_file)
        print('Results:', '         \n'.join('%s_%s.html' % (p, time_stamp) for p in profiles))
        print('\nConformance test', 'passed.' if all_passed else 'failed!', '\n')
        sys.exit(0 if all_passed else 1)
    except KeyboardInterrupt:
        print('\rExiting. Please wait while workers shut down ...')
    finally:
        if not args.keep_files:
            removedir(work_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LC3plus conformance tool - checks if your version of the LC3plus cod'
                                                 'ec is conform to the reference provided by Fraunhofer & Ericsson.')
    parser.add_argument('-keep', action='store_true', dest='keep_files', help="Don't delete workdir at end of test")
    parser.add_argument('-system_sox', action='store_true', help='Use system sox')
    parser.add_argument('-v', action='store_true', dest='verbose', help='Activate verbose output')
    parser.add_argument('-w', dest='workers', type=int, default=os.cpu_count(), help='Number of worker threads')
    parser.add_argument('config', help='Conformance config file')
    args = parser.parse_args()

    if args.system_sox:
        SOX_EXE = 'sox'

    main(args)
