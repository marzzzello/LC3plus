/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

// Provide stuff that's missing on non-linux systems

#ifndef SHIMS_H
#define SHIMS_H

// older visual studio don't have stdint.h and stdbool.h
#if defined _MSC_VER && _MSC_VER < 1800
#error "Visual Studio 2013 or later required"
#endif

// commonly used stdlib headers
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// linux
#if defined __linux__ || defined __CYGWIN__

#include <pthread.h>
#include <endian.h>


// osx
#elif defined __APPLE__

#include <pthread.h>
#include <libkern/OSByteOrder.h>

#define le16toh(x) OSSwapHostToLittleInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)


// windows
#elif defined _WIN32

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

// vs compiler specific fixes
#if defined _MSC_VER

// inline keyword
#define inline __inline

// restrict keyword
#if _MSC_VER >= 1900
#define restrict __restrict
#else
#define restrict
#endif

// disable truncation from 'double' to 'float' warnings
#pragma warning(disable : 4305 4244)

#endif // _MSC_VER

// endian - assuming little endian architecture
#define le16toh(x) ((uint16_t)(x))
#define htole16(x) ((uint16_t)(x))
#define le32toh(x) ((uint32_t)(x))
#define htole32(x) ((uint32_t)(x))

// strcasecmp
#define strcasecmp(s1, s2) _stricmp((s1), (s2))

// sloppy pthread simulation
typedef HANDLE pthread_t;

static inline int pthread_create(pthread_t* thread, const void* attr, void* (*start_routine)(void*), void* arg)
{
    assert(attr == NULL);
    *thread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)start_routine, arg, 0, NULL);
    return *thread ? 0 : -1;
}

static inline int pthread_join(pthread_t thread, void** retval)
{
    assert(retval == NULL);
    DWORD err = WaitForSingleObject(thread, INFINITE);
    CloseHandle(thread);
    return err ? -1 : 0;
}

#else
#error "This system doesn't appear to be supported"
#endif

#endif
