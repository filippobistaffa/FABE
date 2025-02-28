/*
 * memory.c: safer memory allocation
 *
 * Copyright (C) 2008-2016 Daniel P. Berrange
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 */


#ifndef MEMORY_H_
#define MEMORY_H_

#include <config.h>
#include <stdlib.h>
#include <stddef.h>

#include "internal.h"

/**
 * ALLOC:
 * @ptr: pointer to hold address of allocated memory
 *
 * Allocate sizeof(*ptr) bytes of memory and store
 * the address of allocated memory in 'ptr'. Fill the
 * newly allocated memory with zeros.
 *
 * Returns -1 on failure, 0 on success
 */
#define ALLOC(ptr) mem_alloc_n(&(ptr), sizeof(*(ptr)), 1)

/**
 * ALLOC_N:
 * @ptr: pointer to hold address of allocated memory
 * @count: number of elements to allocate
 *
 * Allocate an array of 'count' elements, each sizeof(*ptr)
 * bytes long and store the address of allocated memory in
 * 'ptr'. Fill the newly allocated memory with zeros.
 *
 * Returns -1 on failure, 0 on success
 */
#define ALLOC_N(ptr, count) mem_alloc_n(&(ptr), sizeof(*(ptr)), (count))

/**
 * REALLOC_N:
 * @ptr: pointer to hold address of allocated memory
 * @count: number of elements to allocate
 *
 * Re-allocate an array of 'count' elements, each sizeof(*ptr)
 * bytes long and store the address of allocated memory in
 * 'ptr'. Fill the newly allocated memory with zeros
 *
 * Returns -1 on failure, 0 on success
 */
#define REALLOC_N(ptr, count) mem_realloc_n(&(ptr), sizeof(*(ptr)), (count))

/**
 * FREE:
 * @ptr: pointer holding address to be freed
 *
 * Free the memory stored in 'ptr' and update to point
 * to NULL.
 */
#define FREE(ptr)                               \
  do {                                          \
    free(ptr);                                  \
    (ptr) = NULL;                               \
  } while(0)

/* Return 1 if an array of N objects, each of size S, cannot exist due
   to size arithmetic overflow.  S must be positive and N must be
   nonnegative.  This is a macro, not an inline function, so that it
   works correctly even when SIZE_MAX < N.

   By gnulib convention, SIZE_MAX represents overflow in size
   calculations, so the conservative dividend to use here is
   SIZE_MAX - 1, since SIZE_MAX might represent an overflowed value.
   However, malloc (SIZE_MAX) fails on all known hosts where
   sizeof (ptrdiff_t) <= sizeof (size_t), so do not bother to test for
   exactly-SIZE_MAX allocations on such hosts; this avoids a test and
   branch when S is known to be 1.  */
# define xalloc_oversized(n, s) \
    ((size_t) (sizeof (ptrdiff_t) <= sizeof (size_t) ? -1 : -2) / (s) < (n))

/**
 * mem_alloc_n:
 * @ptrptr: pointer to pointer for address of allocated memory
 * @size: number of bytes to allocate
 * @count: number of elements to allocate
 *
 * Allocate an array of memory 'count' elements long,
 * each with 'size' bytes. Return the address of the
 * allocated memory in 'ptrptr'.  The newly allocated
 * memory is filled with zeros.
 *
 * Returns -1 on failure to allocate, zero on success
 */
__attribute__((always_inline))
static inline int mem_alloc_n(void *ptrptr, size_t size, size_t count)
{
    if (AUGEAS_UNLIKELY(size == 0 || count == 0)) {
        *(void **)ptrptr = NULL;
        return 0;
    }

    *(void**)ptrptr = calloc(count, size);
    if (AUGEAS_UNLIKELY(*(void**)ptrptr == NULL))
        return -1;
    return 0;
}

/**
 * virReallocN:
 * @ptrptr: pointer to pointer for address of allocated memory
 * @size: number of bytes to allocate
 * @count: number of elements in array
 *
 * Resize the block of memory in 'ptrptr' to be an array of
 * 'count' elements, each 'size' bytes in length. Update 'ptrptr'
 * with the address of the newly allocated memory. On failure,
 * 'ptrptr' is not changed and still points to the original memory
 * block. The newly allocated memory is filled with zeros.
 *
 * Returns -1 on failure to allocate, zero on success
 */
__attribute__((always_inline))
static inline int mem_realloc_n(void *ptrptr, size_t size, size_t count)
{
    void *tmp;
    if (AUGEAS_UNLIKELY(size == 0 || count == 0)) {
        free(*(void **)ptrptr);
        *(void **)ptrptr = NULL;
        return 0;
    }
    if (AUGEAS_UNLIKELY(xalloc_oversized(count, size))) {
        errno = ENOMEM;
        return -1;
    }
    tmp = realloc(*(void**)ptrptr, size * count);
    if (AUGEAS_UNLIKELY(!tmp))
        return -1;
    *(void**)ptrptr = tmp;
    return 0;
}

#endif /* __VIR_MEMORY_H_ */
