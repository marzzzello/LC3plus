/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


static int Dec_LC3_Channel(LC3_Dec *decoder, int channel, int bits_per_sample, UWord8 *bs_in, void *s_out, Word16 bfi,
                           Word8 *scratchBuffer)
{
    Word16 scale;
    Word32 offset;
    Word16 fill_bits;
    Word16 nf_seed, gg_idx, fac_ns_idx, q_fx_exp = 0;
    Word16 bp_side, mask_side;
    Word16 tns_numfilters, lsbMode, lastnz, BW_cutoff_idx, BW_cutoff_idx_nf;
    Word16 zero_frame;
#ifdef ENABLE_RFRAME
    Word16 rframe = 0;
#endif
    Word16 ltpf_idx[3];
    Word16 spec_inv_idx = 0;
    Counter i;

    /* Buffers */
    Word16 *  int_scf_fx_exp, tns_order[TNS_NUMFILTERS_MAX];
    Word16 *  resBitBuf;
    Word16 *  sqQdec, *int_scf_fx, *x_fx, *indexes, *scf_q;
    Word32 *  L_scf_idx;
    Word32 *  q_d_fx;
    Word8 *   currentScratch;
    DecSetup *h_DecSetup = decoder->channel_setup[channel];

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Counter i;
        Word16  scale;
        Word32  offset;
        Word16  fill_bits;
        Word16  nf_seed, gg_idx, fac_ns_idx, q_fx_exp;
        Word16  bp_side, mask_side;
        Word16  tns_numfilters, lsbMode, lastnz, BW_cutoff_idx, BW_cutoff_idx_nf;
        Word16  zero_frame;
        Word16  ltpf_idx[3];
#ifdef ENABLE_RFRAME
        Word16 rframe;
#endif
        Word16 spec_inv_idx;

        /* Buffers */
        Word16 *int_scf_fx_exp, tns_order[TNS_NUMFILTERS_MAX];
        Word16 *resBitBuf;
        Word16 *sqQdec, *int_scf_fx, *x_fx, *indexes, *scf_q;
        Word32 *L_scf_idx;
        Word32 *q_d_fx;
        Word8 * currentScratch;
    };
    Dyn_Mem_In("Dec_LC3_Channel", sizeof(struct _dynmem));
#endif

#ifdef DISABLE_PLC
    UNUSED(decoder->plcMeth);
#endif

    /* BUFFER INITIALISATION. Some buffers may overlap since they are not used in the whole decoding process */
    q_d_fx = scratchAlign(scratchBuffer, 0); /* Size = 4 * MAX_LEN bytes */
    resBitBuf =
        scratchAlign(q_d_fx, sizeof(*q_d_fx) * decoder->frame_length); /* Size = 2 * NPRM_RESQ = 2 * MAX_LEN bytes */
    indexes = scratchAlign(
        resBitBuf, sizeof(*resBitBuf) * decoder->frame_length); /* Size = 2 * TNS_NUMFILTERS_MAX * MAXLAG = 32 bytes */
    L_scf_idx      = scratchAlign(indexes, sizeof(*indexes) * TNS_NUMFILTERS_MAX *
                                          MAXLAG); /* Size = 4 * SCF_MAX_PARAM = 28 bytes -> aligned to 32 bytes */
    sqQdec         = scratchAlign(L_scf_idx, sizeof(*L_scf_idx) * (SCF_MAX_PARAM));   /* Size = 2 * MAX_LEN bytes */
    scf_q          = scratchAlign(sqQdec, sizeof(*sqQdec) * (decoder->frame_length)); /* Size = 2 * M = 32 bytes */
    int_scf_fx_exp = scratchAlign(scf_q, sizeof(*scf_q) * M); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    int_scf_fx     = scratchAlign(int_scf_fx_exp,
                              sizeof(*int_scf_fx_exp) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    currentScratch = scratchAlign(int_scf_fx, sizeof(*int_scf_fx) * MAX_BANDS_NUMBER); /* Size = 4 * MAX_LEN */
    x_fx =
        scratchAlign(q_d_fx, sizeof(*q_d_fx) * decoder->frame_length); /* Size = 2 * (MAX_LEN + MDCT_MEM_LEN_MAX) = 2
                                                                        * MAX_LEN + 1.25 * MAX_LEN = 3.25 * MAX_LEN */

#ifdef DISABLE_PLC
    memset(q_d_fx, 0, decoder->frame_length * sizeof(*q_d_fx));
#endif

    BASOP_sub_start("Decoder");

#ifdef ENABLE_RFRAME
    IF (sub(bfi, 3) == 0)
    {
        bfi = 2;  move16();
        rframe = 1;  move16();
    }
#endif

   if ( bfi != 1) 
   {
      BASOP_sub_sub_start("Dec(bfi=0)");
   }
   else
   {
      BASOP_sub_sub_start("Dec(bfi=1)");
   }

    BASOP_sub_start("Entropy dec");
    IF (sub(bfi, 1) != 0)
    {
        processDecoderEntropy_fx(bs_in, &bp_side, &mask_side, h_DecSetup->total_bits, decoder->yLen,
                                 decoder->fs_idx, decoder->BW_cutoff_bits, &tns_numfilters, &lsbMode, &lastnz, &bfi,
                                 tns_order, &fac_ns_idx, &gg_idx, &BW_cutoff_idx, ltpf_idx, L_scf_idx,
                                 decoder->frame_dms);
        BW_cutoff_idx_nf = BW_cutoff_idx;  move16();
    }
    BASOP_sub_end(); /* Entropy dec */

    BASOP_sub_start("Ari dec");
    IF (sub(bfi, 1) != 0)
    {
        processAriDecoder_fx(bs_in, &bp_side, &mask_side, h_DecSetup->total_bits, decoder->yLen, decoder->fs_idx,
                             h_DecSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &bfi, tns_order,
                             fac_ns_idx, gg_idx, decoder->frame_dms,
                             decoder->n_pc, decoder->be_bp_left, decoder->be_bp_right, 0, &spec_inv_idx, &scale,
                             &fill_bits, sqQdec, &nf_seed, resBitBuf, indexes, &zero_frame, currentScratch);
#ifdef ENABLE_RFRAME
        test();
        IF (sub(rframe, 1) == 0 && zero_frame == 0)
        {
            bfi = 2;  move16();
            spec_inv_idx = s_max(lastnz, BW_cutoff_bin_all[BW_cutoff_idx]);  move16();
        }
#endif
        IF (bfi == 0)
        {
            processAriDecoderScaling_fx(sqQdec, decoder->yLen, q_d_fx, &q_fx_exp);
        }
    }
    BASOP_sub_end(); /* Ari dec */

#ifdef BE_MOVED_STAB_FAC
    BASOP_sub_start("SnsQuantScfDec");
    IF (sub(bfi, 1) != 0)
    {
        /* currentScratch Size = 96 bytes */
        processSnsQuantizeScfDecoder_fx(L_scf_idx, scf_q, currentScratch);
    }
    BASOP_sub_end();

    BASOP_sub_start("PLC::ComputeStabFac");
    if (h_DecSetup->plcAd)
    {
        processPLCcomputeStabFac_main(scf_q, h_DecSetup->plcAd->old_scf_q, h_DecSetup->plcAd->old_old_scf_q,
                                      bfi, h_DecSetup->prev_bfi, h_DecSetup->prev_prev_bfi, &h_DecSetup->plcAd->stab_fac);
    }
    BASOP_sub_end();
#endif

    BASOP_sub_start("Partial Concealment");
    IF (sub(bfi, 1) != 0)
    {
        scale = 32767; move16();

        IF (h_DecSetup->plcAd)
        {
            scale = h_DecSetup->plcAd->stab_fac;
        }

        processPCmain_fx(rframe, &bfi, h_DecSetup->prev_bfi, decoder->yLen, decoder->frame_dms,
                         h_DecSetup->q_old_res_fx, &h_DecSetup->q_old_res_fx_exp, sqQdec,
                         h_DecSetup->q_old_d_fx, spec_inv_idx, ltpf_idx[0], scale, q_d_fx, &q_fx_exp,
                         gg_idx, h_DecSetup->quantizedGainOff, &h_DecSetup->prev_gg, &h_DecSetup->prev_gg_e,
                         &BW_cutoff_idx_nf, &h_DecSetup->prev_BW_cutoff_idx_nf, fac_ns_idx, &h_DecSetup->prev_fac_ns_fx,
                         &h_DecSetup->pc_nbLostFramesInRow);
    }
    BASOP_sub_end();

    IF (sub(bfi, 1) != 0)
    {
        BASOP_sub_start("Residual dec");
        processResidualDecoding_fx(q_d_fx, q_fx_exp, decoder->yLen, resBitBuf, fill_bits);
        BASOP_sub_end();

        BASOP_sub_start("Noisefill");
        /* currentScratch Size = 2 * MAX_LEN bytes */
        IF (zero_frame == 0)
        {
            processNoiseFilling_fx(q_d_fx, nf_seed, q_fx_exp, fac_ns_idx, BW_cutoff_idx_nf, decoder->frame_dms,
                                   h_DecSetup->prev_fac_ns_fx, spec_inv_idx, currentScratch);
        }
        BASOP_sub_end();

        BASOP_sub_start("applyGlobalGain");
        processApplyGlobalGain_fx(q_d_fx, &q_fx_exp, decoder->yLen, gg_idx, h_DecSetup->quantizedGainOff);
        BASOP_sub_end();

        BASOP_sub_start("Tns_dec");
        /* currentScratch Size = 48 bytes */
        processTnsDecoder_fx(indexes, q_d_fx, decoder->yLen, tns_order, &q_fx_exp, BW_cutoff_idx, decoder->frame_dms,
                             currentScratch);
        BASOP_sub_end();

#ifndef BE_MOVED_STAB_FAC
        BASOP_sub_start("SnsQuantScfDec");
        /* currentScratch Size = 96 bytes */
        processSnsQuantizeScfDecoder_fx(L_scf_idx, scf_q, currentScratch);
        BASOP_sub_end();
#endif

        BASOP_sub_start("SnsInterpScfDec");
        /* currentScratch Size = 128 bytes */
        processSnsInterpolateScf_fx(scf_q, int_scf_fx, int_scf_fx_exp, 0, decoder->bands_number, currentScratch);
        BASOP_sub_end();

        BASOP_sub_start("Mdct shaping_dec");
        processScfScaling(int_scf_fx_exp, decoder->bands_number, &q_fx_exp);
        processMdctShaping_fx(q_d_fx, int_scf_fx, int_scf_fx_exp, decoder->bands_offset, decoder->bands_number);
        BASOP_sub_end();
        /* end int_scf_fx */
    }

    BASOP_sub_start("PLC::Main");
    /* currentScratch Size = 2 * MAX_LGW + 8 * MAX_LPROT + 12 * MAX_L_FRAME */
    processPLCmain_fx(decoder->plcMeth, &h_DecSetup->concealMethod, &h_DecSetup->nbLostFramesInRow, bfi,
                      h_DecSetup->prev_bfi, decoder->frame_length, decoder->la_zeroes, decoder->W_fx, x_fx,
                      h_DecSetup->stDec_ola_mem_fx, &h_DecSetup->stDec_ola_mem_fx_exp, h_DecSetup->q_old_d_fx,
                      &h_DecSetup->q_old_fx_exp, q_d_fx, &q_fx_exp, decoder->yLen, decoder->fs_idx,
                      decoder->bands_offset, &h_DecSetup->plc_damping, h_DecSetup->ltpf_mem_pitch_int,
                      h_DecSetup->ltpf_mem_pitch_fr, &h_DecSetup->ns_cum_alpha, &h_DecSetup->ns_seed, h_DecSetup->plcAd,
                      decoder->frame_dms, currentScratch);
    BASOP_sub_end();

#ifdef NONBE_PLC4_ADAP_DAMP
    BASOP_sub_start("PLC/PC::DampingScrambling");
    if (h_DecSetup->plcAd)
    {
        processPLCDampingScrambling_main_fx(bfi, h_DecSetup->concealMethod, h_DecSetup->nbLostFramesInRow,
                                            h_DecSetup->pc_nbLostFramesInRow, &h_DecSetup->ns_seed, &h_DecSetup->pc_seed,
                                            h_DecSetup->ltpf_mem_pitch_int, ltpf_idx[0], q_d_fx, &q_fx_exp, h_DecSetup->q_old_d_fx,
                                            &h_DecSetup->q_old_fx_exp, decoder->yLen, h_DecSetup->plcAd->stab_fac, decoder->frame_dms,
                                            &h_DecSetup->plcAd->cum_fading_slow, &h_DecSetup->plcAd->cum_fading_fast,
                                            &h_DecSetup->plc_damping, spec_inv_idx);
    }
    BASOP_sub_end();
#endif

    BASOP_sub_start("Imdct");
    /* currentScratch Size = 4 * MAX_LEN */
    ProcessingIMDCT(q_d_fx, &q_fx_exp, decoder->W_fx, h_DecSetup->stDec_ola_mem_fx, &h_DecSetup->stDec_ola_mem_fx_exp,
                    x_fx, decoder->W_size, decoder->frame_length, decoder->stDec_ola_mem_fx_len, decoder->frame_dms,
                    h_DecSetup->concealMethod, bfi, h_DecSetup->prev_bfi, h_DecSetup->nbLostFramesInRow,
                    h_DecSetup->plcAd,
                    currentScratch);
    BASOP_sub_end();

   
       BASOP_sub_start("PLC::Update");
    
  

    processPLCupdate_fx(h_DecSetup->plcAd, x_fx, q_fx_exp, h_DecSetup->concealMethod, decoder->frame_length,
                        decoder->fs_idx, &h_DecSetup->nbLostFramesInRow, &h_DecSetup->prev_prev_bfi, &h_DecSetup->prev_bfi,
                        bfi, scf_q, h_DecSetup->stDec_ola_mem_fx, h_DecSetup->stDec_ola_mem_fx_exp, &h_DecSetup->ns_cum_alpha);
    BASOP_sub_end();

#ifdef LTPF_DISABLE_FILTERING   
    ltpf_idx[0] = 0;
    ltpf_idx[1] = 0;
    ltpf_idx[2] = 0;
    h_DecSetup->ltpf_mem_active=0;
#endif


    BASOP_sub_start("LtpfDec");
    /* currentScratch Size = 0.5 * MAX_LEN + 20 bytes */
    process_ltpf_decoder_fx(&q_fx_exp, decoder->frame_length, decoder->ltpf_mem_x_len, decoder->fs_idx,
                            decoder->ltpf_mem_y_len, &h_DecSetup->ltpf_mem_e, x_fx, h_DecSetup->ltpf_mem_x, x_fx,
                            h_DecSetup->ltpf_mem_y, ltpf_idx[0], ltpf_idx[1], ltpf_idx[2],
                            &h_DecSetup->ltpf_mem_pitch_int, &h_DecSetup->ltpf_mem_pitch_fr, &h_DecSetup->ltpf_mem_gain,
                            &h_DecSetup->ltpf_mem_active, h_DecSetup->ltpf_scale_fac_idx, bfi,
                            h_DecSetup->concealMethod,
                            h_DecSetup->plc_damping, &h_DecSetup->ltpf_mem_scale_fac_idx, currentScratch);
    BASOP_sub_end();

    BASOP_sub_start("Output scaling");
    {
        scale  = sub(sub(31 + 16, bits_per_sample), q_fx_exp);
        offset = L_shr_sat(32768, sub(16, scale));
        IF (bits_per_sample == 16)
        {
            scale = sub(15, q_fx_exp);
            FOR (i = 0; i < decoder->frame_length; i++)
            {
                ((Word16 *)s_out)[i] = round_fx_sat(L_shr_sat(L_deposit_h(x_fx[i]), scale)); move16();
            }
        }
        ELSE
        {
            FOR (i = 0; i < decoder->frame_length; i++)
            {
                ((Word32 *)s_out)[i] = L_shr_sat(L_add_sat(L_deposit_h(x_fx[i]), offset), scale); move32();
            }
        }
    }
    BASOP_sub_end(); /* Output scaling */
    
    BASOP_sub_sub_end();

    BASOP_sub_end(); /* Decoder */


#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    return bfi;
}

/* num_bytes = 0 -> bad frame */
LC3_Error Dec_LC3(LC3_Dec *decoder, UWord8 *input, int num_bytes, void **output, int bits_per_sample, void *scratch, int bfi_ext)
{
    int       ch = 0, bfi = bfi_ext;
    LC3_Error err = LC3_OK;
    int       fec_num_bytes;
    int       lc3_num_bytes;
    int       lc3_channel_num_bytes;
    int       channel_bfi, out_bfi;
    Word16    channel_epmr;

    if (bfi == 0)
    {
        bfi = !num_bytes;
    }

    if (decoder->ep_enabled)
    {
        decoder->combined_channel_coding = decoder->channels > 1 && num_bytes <= 160;

        if (decoder->combined_channel_coding)
        {
            fec_num_bytes = num_bytes;

            BASOP_sub_start("fec_dec");

            decoder->error_report =
                fec_decoder(input, fec_num_bytes, &lc3_num_bytes, &decoder->epmr, decoder->combined_channel_coding,
                            &decoder->n_pccw, &bfi, &decoder->be_bp_left, &decoder->be_bp_right, &decoder->n_pc, &decoder->m_fec, scratch);

            BASOP_sub_end();

            for (ch = 0; ch < decoder->channels; ch++)
            {
                lc3_channel_num_bytes = lc3_num_bytes / decoder->channels + (ch < (lc3_num_bytes % decoder->channels));


                if (bfi != 1 && lc3_channel_num_bytes != decoder->channel_setup[ch]->last_size)
                {
                    err = update_dec_bitrate(decoder, ch, lc3_channel_num_bytes);

                    if (err)
                        return err;

                    decoder->channel_setup[ch]->last_size = lc3_channel_num_bytes;
                }

                bfi = Dec_LC3_Channel(decoder, ch, bits_per_sample, input, output[ch], bfi, scratch);
                input += decoder->channel_setup[ch]->targetBytes;
            }
        }
        else
        {
            decoder->epmr = 12;
            out_bfi      = 0;

            for (ch = 0; ch < decoder->channels; ch++)
            {
                fec_num_bytes = num_bytes / decoder->channels + (ch < (num_bytes % decoder->channels));

                BASOP_sub_start("fec_dec");

                channel_bfi = bfi;

                decoder->error_report = fec_decoder(input, fec_num_bytes, &lc3_num_bytes, &channel_epmr,
                                                    decoder->combined_channel_coding, &decoder->n_pccw, &channel_bfi,
                                                    &decoder->be_bp_left, &decoder->be_bp_right, &decoder->n_pc, &decoder->m_fec, scratch);

                BASOP_sub_end();

                decoder->epmr = MIN(decoder->epmr, channel_epmr);


#ifdef ENABLE_PADDING
				if (channel_bfi != 1)
                {
                    Word16 padding_len, np_zero;

                    if (paddingDec_fx(input, shl(lc3_num_bytes, 3), decoder->yLen, decoder->BW_cutoff_bits, decoder->ep_enabled, &padding_len, &np_zero))
                    {
                    	channel_bfi = 1;
                    }

                    input         = input + np_zero;
                    decoder->n_pc = s_max(decoder->n_pc - (2 * np_zero), 0);

                    if (channel_bfi == 2)
                    {
                        if (decoder->be_bp_right < (8*np_zero))
                        {
                            channel_bfi = 0;
                            decoder->be_bp_left = -1;
                            decoder->be_bp_right = -1;
						}
                        else
                        {
                            decoder->be_bp_right = decoder->be_bp_right - (8 * np_zero);
                            decoder->be_bp_left  = s_max(decoder->be_bp_left - (8 * np_zero), 0);
						}
					}

                    lc3_num_bytes = lc3_num_bytes - padding_len;
                }
#endif

                if (channel_bfi != 1 && lc3_num_bytes != decoder->channel_setup[ch]->last_size)
                {
                    err = update_dec_bitrate(decoder, ch, lc3_num_bytes);
                    if (err)
                        return err;

                    decoder->channel_setup[ch]->last_size = lc3_num_bytes;
                }

                channel_bfi = Dec_LC3_Channel(decoder, ch, bits_per_sample, input, output[ch], channel_bfi, scratch);

                out_bfi |= channel_bfi;
                input += fec_num_bytes;
            }

            bfi = out_bfi & 1;
        }
    }
    else
    {
        for (ch = 0; ch < decoder->channels; ch++)
        {
            lc3_num_bytes = num_bytes / decoder->channels + (ch < (num_bytes % decoder->channels));

#ifdef ENABLE_PADDING
            if (bfi != 1)
            {
                Word16 padding_len, np_zero;

                if (paddingDec_fx(input, shl(lc3_num_bytes, 3), decoder->yLen, decoder->BW_cutoff_bits, decoder->ep_enabled, &padding_len, &np_zero)){
                    bfi = 1;
                }

                lc3_num_bytes = lc3_num_bytes - padding_len;
                if (lc3_num_bytes < 20 || lc3_num_bytes > LC3_MAX_BYTES) {
                    bfi = 1;    /* mark frame as broken if frame sizeif below the minimum of 20 bytes */
                }

            }
#endif 
			
            if (bfi != 1 && lc3_num_bytes != decoder->channel_setup[ch]->last_size)
            {
                err = update_dec_bitrate(decoder, ch, lc3_num_bytes);
                if (err)
                    return err;
                decoder->channel_setup[ch]->last_size = lc3_num_bytes;
            }

            bfi = Dec_LC3_Channel(decoder, ch, bits_per_sample, input, output[ch], bfi, scratch);
            input += decoder->channel_setup[ch]->targetBytes;
        }
    }

    return bfi == 1 ? LC3_DECODE_ERROR : LC3_OK;
}
