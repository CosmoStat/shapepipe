
/*----------------------------------------------------------------------------
 * 									E.S.O.
 *----------------------------------------------------------------------------
 * File name 	:	memory.h
 * Author 		:	Nicolas Devillard
 * Created on	:	Nov 07, 1995
 * Hardware		:	Sun Sparc 20
 * Software		:	ANSI C under Solaris Unix
 *					Part of ECLIPSE library for Adonis
 * Description	:	memory handling and swapping routines
 *--------------------------------------------------------------------------*/
 
#ifndef _VMS_H_
#define _VMS_H_

#include <stdlib.h>


char *mem_alloc_buffer(size_t  Nelem);
void mem_free_buffer(char *Ptr);





#ifdef LARGE_BUFF

#define DEF_VMS_SIZE 4  /* megabytes */
#define DEF_VMS_DIR "."

#define EXPORT_MEMORY_H

/*----------------------------------------------------------------------------
 * 								Includes
 *--------------------------------------------------------------------------*/

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string>
#ifndef MACOS
#include <malloc.h>
#endif
#include <signal.h>
#include <cerrno>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <sys/resource.h>

#define e_error printf
#define e_warning printf

#ifdef  __cplusplus
extern "C" {
#endif


typedef struct _MMAP_FILE_ {
	int             fd ;
	void  *         buf ;
	unsigned long   size ;
	char            filename[1024] ;
} mmap_file ;


/*----------------------------------------------------------------------------
 * 								Defines	
 *--------------------------------------------------------------------------*/


#define ONE_MEG     (size_t)1048576

#ifndef NAME_SIZE
#define NAME_SIZE   512
#endif


/*
 * Allow here extended memory.
 */

#define ENABLE_EXTENDED_MEMORY



/*----------------------------------------------------------------------------
 * 								Macros	
 *--------------------------------------------------------------------------*/

/*
 * Probably the most useful macro: NONNULL(ptr,retval) will exit the
 * current function and return [retval] if ptr is found to be 0
 */

#define NONNULL(ptr, retval) { if (!ptr) return retval ; }

#define ext_malloc(s) 		ext_malloc_x(s,__FILE__,__LINE__)
#define ext_calloc(n,s)		ext_calloc_x(n,s,__FILE__,__LINE__)
#define ext_free(p) 		ext_free_x(p,__FILE__,__LINE__)
#define swap_malloc(s)		swap_malloc_x(s,__FILE__,__LINE__)
#define swap_free(p)		swap_free_x(p,__FILE__,__LINE__)

#define malloc(s) ext_malloc_x(s,__FILE__,__LINE__)
#define calloc(n,s) ext_calloc_x(n,s,__FILE__,__LINE__)
#define free(p) ext_free_x(p,__FILE__,__LINE__)
#define strdup(s) ext_strdup_x(s,__FILE__,__LINE__)

/*----------------------------------------------------------------------------
 *						Function ANSI C prototypes
 *--------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 * Function :   ext_malloc_x()
 * In       :   size_t
 * Out      :   pointer to allocated zone
 * Job      :   allocate memory in RAM or SWAP
 * Notice   :
 *--------------------------------------------------------------------------*/

void *
ext_malloc_x(
    size_t      size,
    char    *   filename,
    int         lineno
) ;


/*----------------------------------------------------------------------------
 * Function :   ext_calloc_x()
 * In       :   size_t
 * Out      :   pointer to allocated (and zero-filled) zone
 * Job      :   allocate memory in RAM or SWAP
 * Notice   :
 *--------------------------------------------------------------------------*/

void *
ext_calloc_x(
	size_t		n_elem,
    size_t      size,
    char    *   filename,
    int         lineno
) ;


/*---------------------------------------------------------------------------
 * Function :   ram_malloc()
 * In       :   Number of Bytes
 * Out      :   pointer to allocated zone
 * Job      :   Allocates memory (RAM)
 * Notice   :   returns NULL in case of error
 *--------------------------------------------------------------------------*/

void *
ram_malloc(
    size_t      size,
    char    *   filename,
    int         lineno
) ;


/*---------------------------------------------------------------------------
 * Function :   ram_calloc()
 * In       :   number of elements, size of each element
 * Out      :   pointer to allocated zone
 * Job      :   Allocates memory (RAM)
 * Notice   :   returns NULL in case of error
 *--------------------------------------------------------------------------*/

void *
ram_calloc(
    size_t      n_elem,
    size_t      size,
    char    *   filename,
    int         lineno
) ;



/*---------------------------------------------------------------------------
 * Function :   ext_free_x()
 * In       :   Pointer to de-allocate
 * Out      :   void
 * Job      :   handles pointer free
 * Notice   :
 *--------------------------------------------------------------------------*/

void
ext_free_x(
    void    *   pointer,
    char    *   filename,
    int         lineno
) ;



/*---------------------------------------------------------------------------
   Function :   swap_malloc_x()
   In       :   Number of bytes to allocate in swap area
   Out      :   pointer to allocated zone
   Job      :   Allocates memory in swap space, returns NULL if cannot
   Notice   :
                memory should be freed with swap_free_monit()
                this function should not be called directly but through
                the swap_malloc() macro
 ---------------------------------------------------------------------------*/

void *
swap_malloc_x(
    size_t      size,
    char    *   filename,
    int         lineno
) ;



/*---------------------------------------------------------------------------
 * Function :   swap_free_x()
 * In       :   SWAP pointer to de-allocate
 * Out      :   void
 * Job      :   handles swap pointer free
 * Notice   :   don't even try to understand...
 *--------------------------------------------------------------------------*/

void
swap_free_x(
    void    *   pointer,
    char    *   filename,
    int         lineno
) ;


/*---------------------------------------------------------------------------
   Function :   ext_strdup_x()
   In       :   char *
   Out      :   char *
   Job      :   copy a string to a newly allocated string
   Notice   :   returned string must be freed by free()
 ---------------------------------------------------------------------------*/

char * ext_strdup_x(char * s, char * file, int lineno) ;



/*----------------------------------------------------------------------------
 * Function :   get_memory_parameter()
 * In       :   allocated character string
 * Out      :   size_t
 * Job      :   get info about memory handling configuration
 * Notice   :   The I/O mapping is the following:

    Input char string               output value

    "max_ram"                       maximum RAM to allocate
    "max_swap"                      maximum SWAP to allocate
    "vm_page_size"                  page size for virtual memory
    "total_alloc"                   total allocated memory so far
    "total_alloc_ram"               total allocated RAM so far
    "total_alloc_swap"              total allocated SWAP so far

 *--------------------------------------------------------------------------*/

size_t
get_memory_parameter(char * s) ;



/*----------------------------------------------------------------------------
 * Function :   set_memory_parameter()
 * In       :   char string, size_t
 * Out      :   void
 * Job      :   set info about memory handling configuration
 * Notice   :   The I/O mapping is the following:

    Input char string               output value

    "max_ram"                       maximum RAM to allocate
    "max_swap"                      maximum SWAP to allocate
    "vm_page_size"                  page size for virtual memory

 *--------------------------------------------------------------------------*/

void
set_memory_parameter(char * s, size_t val) ;



/*---------------------------------------------------------------------------
   Function :   print_memory_parameters()
   In       :   void
   Out      :   void, messages to stderr
   Job      :   print out current memory configuration
   Notice   :
 ---------------------------------------------------------------------------*/

void print_memory_parameters(void) ;



/*---------------------------------------------------------------------------
 * Function :   get_tempfile_name()
 * In       :   void
 * Out      :   newly allocated character string
 * Job      :   return a different file name at each call
 * Notice   :
 *--------------------------------------------------------------------------*/
char *
get_tempfile_name(void) ;


/*---------------------------------------------------------------------------
   Function :   set_tmpdirname()
   In       :   char *
   Out      :   void
   Job      :   sets the name of the temporary directory for VM
   Notice   :   
 ---------------------------------------------------------------------------*/

void
set_tmpdirname(char * s) ;


/*---------------------------------------------------------------------------
 * Function :   memory_status()
 * In       :   void
 * Out      :   prints out memory status on stderr
 * Notice   :
 *--------------------------------------------------------------------------*/
void
memory_status(void) ;


/*----------------------------------------------------------------------------
 * Function :   cleanup_mess()
 * In       :   void
 * Out      :   void
 * Job      :   clean up temporary files created by the current process
 * Notice   :
 *--------------------------------------------------------------------------*/

void cleanup_mess(void) ;


/*----------------------------------------------------------------------------
 * Function :   signal_catch()
 * In       :   signal to catch, function to call in this case
 * Out      :   int 0 if worked, -1 otherwise
 * Job      :   catch an interrupt and launch the appropriate function
 * Notice   :   directly taken from MIDAS
 *--------------------------------------------------------------------------*/

int signal_catch(int sig, void (*f)()) ;


/*---------------------------------------------------------------------------
   Function :   fatal_signal_handler()
   In       :   void
   Out      :   void
   Job      :   prints out a message to stderr and call cleanup
   Notice   :   called from SIGSEGV, SIGBUS, SIGXCPU, SIGXFSZ
 ---------------------------------------------------------------------------*/

void fatal_signal_handler(void) ;


/*----------------------------------------------------------------------------
 * Function :   eclipse_init()
 * In       :   void
 * Out      :   void
 * Job      :   install the atexit() routine and interrupt_catch() to
 *              take care of removing temporary files in case of sudden
 *              process death.
 * Notice   :   by no means compulsory
 *--------------------------------------------------------------------------*/

void eclipse_init(void) ;


/*---------------------------------------------------------------------------
   Function :   enable_virtual_memory()
   In       :   void
   Out      :   int 0 if Ok, anything else otherwise
   Job      :   enables the use of virtual memory, i.e. sets the size of
                one virtual memory page, and checks the available disk
                space to see if max_swap fits.
   Notice   :   should always be called before any use of virtual
                memory.
 ---------------------------------------------------------------------------*/

int enable_virtual_memory(void) ;


/*---------------------------------------------------------------------------
   Function :   mmap_open()
   In       :   filename
   Out      :   pointer to allocated mmap_file structure
   Job      :   mmap the totality of a given file
   Notice   :   should only be closed by mmap_close()
 ---------------------------------------------------------------------------*/


mmap_file *
mmap_open(char * filename) ;


/*---------------------------------------------------------------------------
   Function :   mmap_close()
   In       :   allocated mmapped_file struct pointer
   Out      :   void
   Job      :   munmap and close the file
   Notice   :
 ---------------------------------------------------------------------------*/

void mmap_close(mmap_file * mm) ;


/*---------------------------------------------------------------------------
   Function :   dump_stack()
   In       :   void
   Out      :   void, messages on stderr
   Job      :   dump the current stack to stderr, using dbx
   Notice   :   Probably Solaris specific, with little efforts portable
                other Unixes. From the comp.unix.programmer FAQ.
 ---------------------------------------------------------------------------*/
 
void dump_stack(void);

#ifdef EXPORT_MEMORY_H
size_t filesize(char * filename) ;
#endif

#ifdef  __cplusplus
}
#endif
#endif
#endif
/*--------------------------------------------------------------------------*/
