
/*----------------------------------------------------------------------------
 * 									E.S.O.
 *----------------------------------------------------------------------------
 * File name 	:	VMS.cv
 * Author 	:	Nicolas Devillard
 * Created on	:	Nov 07, 1995
 * Language	:	ANSI C
 *			Virtual Memory Software
 * Description	:	memory handling and swapping routines

	The following module offers functionalities for efficient memory
	handling. The basic goal is to move the upper limit of usable RAM on
	the machine to the amount of available disk space (within limits
	fixed by the OS).
	With this library, it is possible to:

	* log all calls to malloc/calloc/free, logging data such as the
	  allocated amount of space, the calling file and line number, the
	  returned pointer value and the exact data location.
	* check calls consistency: warn when a pointer is freed several
	  times, or when a block has not been deallocated at the end of the
	  program.
	* print out current process memory usage.
	* allows the use of virtual memory, in an invisible way.
	* checks out that temporary files are deleted.
 
        * if LARGE_BUFF compilation option is not define, VMS is not
	  used, and standard free and malloc are user
 *--------------------------------------------------------------------------*/
#ifndef _INC_STDLIB_H
#include <stdlib.h>
#define _INC_STDLIB_H
#endif

#ifndef _IM_GLOB_H_
#include "GlobalInc.h"
#endif
#ifndef _VMS_H_
#include "VMS.h"
#endif
#ifndef _TEMPMEMORY_H
#include "TempMemory.h"
#endif

#ifdef LARGE_BUFF
/* #define the following to get zillions of debug messages */
// #define DEBUG_ALLOC 

#define	ALLOC_CELL_NOALLOC	-1
#define	ALLOC_CELL_RAM		 0
#define ALLOC_CELL_SWAP		 1

#define SRCFILENAME_SZ		256
#define TMPFILENAME_SZ		256

//size_t emem_ramlimit = DEF_VMS_SIZE; 	/* megabytes */
int emem_ramlimit = DEF_VMS_SIZE; 	/* megabytes */
char  emem_tmpdirname[1024] =  DEF_VMS_DIR;

#endif


#ifdef LARGE_BUFF
class VMS_Status {
     public:
       Bool Verbose;
       void set_size(int UserSize=0);
       void set_dir(char * UserName=NULL);
       void print();
       VMS_Status()
       {
          Verbose = False;
          set_size();
	  set_dir();
       }
       ~VMS_Status() {Verbose=False;}
};

void VMS_Status::set_size(int UserSize)
{
   extern int emem_ramlimit;
   if (UserSize < 1)
   {
       int VMSSize=0;
       char *PtrSize;
       PtrSize = (char *) getenv(VMS_SIZE);
       if (PtrSize != NULL)
       {
           VMSSize = atoi(PtrSize);
	   if (VMSSize > 0)   emem_ramlimit = VMSSize;
       }
   }
   else emem_ramlimit = UserSize;
}

void VMS_Status::set_dir(char * UserName)  
{
   extern char emem_tmpdirname[1024];
   if ((UserName != NULL) && (strlen(UserName) > 0))
                           strcpy(emem_tmpdirname, UserName); 
   else
   {
      char *PtrName = (char *) getenv(VMS_DIR);
      if (PtrName != NULL)  strcpy(emem_tmpdirname, PtrName);
   }  
}	

void VMS_Status::print()
{
   extern char emem_tmpdirname[1024];
   extern int emem_ramlimit;
   printf("Virtual memory limit size = %d\n", emem_ramlimit);
   printf("Virtual memory directory = %s\n",  emem_tmpdirname);
}
   
VMS_Status VMSStatus;

// ================== END class VMS Stat ========================	
	 
void vms_init(int UserSize, char * UserName, Bool Verbose)
{
   extern Bool UseVMS;
 
   UseVMS = True;
   // temporary file name
   // if a user file is given, take it
   // else look at CEA_VM_DIR environment variable
   // else keep default name
   VMSStatus.Verbose = Verbose;
   VMSStatus.set_dir(UserName);
   
   // minimum size for using the virtual memory
   // if UserSize is given, take it
   // else look at  CEA_VM_SIZE environment variable
   // else keep default value
   VMSStatus.set_size(UserSize);
   if (Verbose == True) VMSStatus.print();
}
#endif

















/****************************************************************************/

char *mem_alloc_buffer(size_t  Nelem)
{
   char *Ptr = (char *) malloc (Nelem);
   return Ptr;
}

/****************************************************************************/

void mem_free_buffer(char *Ptr)
{
   free (Ptr);
}

/****************************************************************************/




#ifdef LARGE_BUFF
#ifdef  __cplusplus
extern "C" {
#endif
 
/*----------------------------------------------------------------------------
 * 								New types
 *--------------------------------------------------------------------------*/

/* This structure holds all allocated pointers (forward linked list)	*/
typedef struct _MALLOC_CELL_ {
	void 				*	pointer;
	size_t					size;
	char					filename[TMPFILENAME_SZ] ;
	int						lineno ;
	int						flavour ;
	char					swapfilename[SRCFILENAME_SZ] ;
	int						swapfd ;
	struct _MALLOC_CELL_*	next ;
	struct _MALLOC_CELL_*	prev ;
} malloc_cell;
		 
/*----------------------------------------------------------------------------
 * 							Global variables
 *--------------------------------------------------------------------------*/

/* This structure holds information about everything allocated so far	*/
static malloc_cell * malloc_table = NULL ;

/*
 * Global variables. Their visibility is limited to this module.
 * Access and modification of their values is done through associated
 * member function calls.
 */
static	size_t	total_allocated_RAM = 0 ;
static	size_t	total_allocated_SWAP = 0 ;
static	int		nswapfiles = 0 ;
static  int		nalloc_ptrs = 0 ;

/* 
 * a register is kept to know how many temporary files have been opened
 * it is only used to generate unique file names
 */
static 	int		file_reg = 0 ;

/*
 * The following parameters set up the limits for this memory module.
 * They are defined as global variables to this module and should not be
 * seen outside except through member functions get/set_memory_parameter().
 */

static size_t emem_vm_limit			= 1024 ;	/* megabytes */
static size_t emem_vm_pagesize		= 1024 ;	/* kilobytes */
//extern size_t emem_ramlimit;
//extern char  emem_tmpdirname[1024];

 /*
main(int argc, char *argv[])
{
   char *P;
   unsigned int N = 0;
   unsigned int M = 0;
   if (argc == 2) N = atoi (argv[1]);
   else 
   {
      printf("USAGE: memory N_Megabits\n");
      exit(0);
   }
   M = 1000000 * N;
   P = (char *) mem_alloc_buffer (M);
   if (P != NULL)
   {
      printf("allocation of %d Mega bytes OK\n", N);
   }
   else
   {
      printf("failed for allocation of %d Mega bytes\n", N);
   }
   mem_free_buffer (P);
   exit(0);
}
*/


#ifdef ENABLE_EXTENDED_MEMORY
#undef malloc
#undef calloc
#undef free
#undef strdup
#endif


/*----------------------------------------------------------------------------
 *							Function codes
 *--------------------------------------------------------------------------*/
#define DF_CMD "df -k"

#ifdef _HPUX_SOURCE
#undef DF_CMD
#define DF_CMD "bdf"
#endif

/*---------------------------------------------------------------------------
   Function	:	test_write_permission()
   In 		:	path
   Out 		:	int: 0 or 1
   Job		:	find out if the current user has write permission on the
				provided path.
   Notice	:
 ---------------------------------------------------------------------------*/


int test_write_permission(char * path)
{
	return !access(path, W_OK);
}


/*---------------------------------------------------------------------------
   Function	:	get_avail_kbytes()
   In 		:	path
   Out 		:	unsigned long
   Job		:	find out how many kbytes are available on the path
   Notice	:	outsourced to 'df' or 'bdf'
 ---------------------------------------------------------------------------*/


unsigned long get_avail_kbytes(char * path)
{
    FILE            *   df_prog ;
    char                line[1024] ;
    char                cmd[1024] ;
    char                filesystem[256] ;
    unsigned long       kbytes ;
    unsigned long       used ;
    unsigned long       avail ;
    int                 capacity ;
    char                df_path[1024] ;

    sprintf(cmd, "%s %s", DF_CMD, path) ;
    df_prog = popen(cmd, "r") ;
    if (df_prog == NULL) {
        printf("errors launching %s: aborting\n", cmd) ;
        return 0 ;
    }

    fgets(line, 1024, df_prog) ;
    fgets(line, 1024, df_prog) ;
	if (pclose(df_prog)==-1) perror("pclose");

    if ((int)strlen(line)<1) {
        e_error("reading output from %s: aborting", cmd) ;
        return 0 ;
    }
    sscanf(line, "%s %ld %ld %ld %d%% %s",
            filesystem,
            &kbytes,
            &used,
            &avail,
            &capacity,
            df_path) ;
    return avail;
}

/*---------------------------------------------------------------------------
 * Function	:	get_tempfile_name()
 * In 		:	void
 * Out 		:	pointer to static character string	
 * Job		:	return a different file name at each call
 * Notice	:	no need to deallocate, it is static memory
 *--------------------------------------------------------------------------*/
char *
get_tempfile_name(void)
{
	static char tmpfilename[TMPFILENAME_SZ] ;
	file_reg++ ;
	(void)sprintf(tmpfilename, "%s/vmswap_%05ld_%05x",
				  emem_tmpdirname,
				  (long)getpid(),
//                                (long)file_reg) ;				  
				  (unsigned int)file_reg) ;
	/* strcpy(tmpfilename, "/dev/zero") ; */
	return tmpfilename ;
}


/*---------------------------------------------------------------------------
   Function	:	set_tmpdirname()
   In 		:	char *	
   Out 		:	void
   Job		:	sets the name of the temporary directory for VM
   Notice	:
 ---------------------------------------------------------------------------*/

void
set_tmpdirname(char * s)
{
	(void)strncpy(emem_tmpdirname, s, TMPFILENAME_SZ-1) ;
	return ;
}




/*----------------------------------------------------------------------------
 * Function	:	ext_malloc_x()
 * In 		:	size_t	
 * Out 		:	pointer to allocated zone
 * Job		:	allocate memory in RAM or SWAP
 * Notice	:
 *--------------------------------------------------------------------------*/

void *ext_malloc_x(
	size_t		size,
	char	*	filename,
	int			lineno)
{
	void	*	pointer ;
          
	/*
	 * See if we need to allocate swap for the requested size.
	 */
 
			    
        if (size > emem_vm_pagesize)
	{
 		if ((size+total_allocated_RAM)>=emem_ramlimit*ONE_MEG) 
		{
        // printf(" size+total_allocated_RAM = %d,  emem_ramlimit*ONE_MEG = %d\n", 
        //                    size+total_allocated_RAM , emem_ramlimit*ONE_MEG);		
			pointer = swap_malloc_x(size, filename, lineno) ;
			return pointer ;
		}
	}

	/*
	 * Otherwise allocate in RAM
	 */
	pointer = ram_malloc(size, filename, lineno) ;
	return pointer ;
}


/*----------------------------------------------------------------------------
 * Function	:	ext_calloc_x()
 * In 		:	size_t	
 * Out 		:	pointer to allocated (and zero-filled) zone
 * Job		:	allocate memory in RAM or SWAP
 * Notice	:
 *--------------------------------------------------------------------------*/

void *
ext_calloc_x(
	size_t		n_elem,
	size_t		size,
	char	*	filename,
	int			lineno
)
{
	size_t		reqsize ;
	void	*	pointer ;

	/*
	 * See if we need to allocate swap for the requested size.
	 */
	reqsize = size * n_elem ;
	if (reqsize  > emem_vm_pagesize) {
		if ((reqsize+total_allocated_RAM)>emem_ramlimit*ONE_MEG) {
			pointer = swap_malloc_x(reqsize, filename, lineno) ;
			return pointer ;
		}
	}
	/*
	 * Otherwise allocate in RAM
	 */
	pointer = ram_calloc(n_elem, size, filename, lineno) ;
	return pointer ;
}



/*---------------------------------------------------------------------------
 * Function	:	ram_malloc()
 * In 		:	Number of Bytes
 * Out 		:	pointer to allocated zone
 * Job		:	Allocates memory (RAM)
 * Notice	:	returns NULL in case of error
 *--------------------------------------------------------------------------*/

void *
ram_malloc(
	size_t		size,
	char	*	filename,
	int			lineno
)
{
	void		*pointer ;
	malloc_cell	*new_cell ;
	malloc_cell	*curr_cell ;

	pointer = (void*)malloc(size) ;
#ifdef DEBUG_ALLOC
	(void)fprintf(stderr,
				  "-> %p = RAM malloc(%ld) in %s (%d)\n",
				  pointer, size, filename, lineno) ;
	(void)fflush(stderr) ;
#endif
	if (pointer == NULL) {
		e_error("failed to allocate RAM %ld bytes in %s (%d)\n", 
				(long)size, filename, lineno) ;
		e_error("total RAM allocated so far : %d bytes\n", total_allocated_RAM);
		e_error("exiting program now\n") ;
		exit(-2000) ;
	}

	total_allocated_RAM += size ;	
	nalloc_ptrs ++ ;
	new_cell = (malloc_cell*)malloc(sizeof(malloc_cell)) ;
	if (new_cell == (malloc_cell*)NULL) {
		e_error("fatal error in memory : life functions terminated\n") ;
		exit(2001) ;
	}
	new_cell->pointer = pointer ;
	new_cell->size = size ;
	strncpy(new_cell->filename, filename, SRCFILENAME_SZ-1) ;
	new_cell->lineno = lineno ; 
	new_cell->flavour = ALLOC_CELL_RAM ; 
	new_cell->next = NULL ;
	new_cell->prev = NULL ;

	/*
	 * Place new memory cell into the structure
	 */
	
	if (malloc_table == NULL) {
		malloc_table = new_cell ;
	} else {
		curr_cell = malloc_table ;
		while (curr_cell->next != NULL) {
			curr_cell = curr_cell->next ;
		}
		curr_cell->next = new_cell ;
		new_cell->prev  = curr_cell ;
	}
	return(pointer) ;
}

/*---------------------------------------------------------------------------
 * Function :   ram_calloc()
 * In       :   number of elements, size of each element
 * Out      :   pointer to allocated zone
 * Job      :   Allocates memory (RAM)
 * Notice   :   returns NULL in case of error
 *--------------------------------------------------------------------------*/

void *
ram_calloc(
	size_t		n_elem,
	size_t		size,
	char	*	filename,
	int			lineno
)
{
    void        *pointer ;
    malloc_cell  *new_cell ;
    malloc_cell  *curr_cell ;

    pointer = (void*)calloc(n_elem, size) ;
#ifdef DEBUG_ALLOC
	(void)fprintf(stderr,
				  "-> %p = RAM calloc(%ld,%ld) in %s (%d)\n",
				  pointer, n_elem, size, filename, lineno) ;
	(void)fflush(stderr) ;
#endif
    if (pointer == NULL) {
        e_error("failed to allocate RAM %ld bytes in %s (%d)\n",
				(long)n_elem * size, filename, lineno) ;
		e_error("total RAM allocated so far: %d bytes\n", total_allocated_RAM);
		e_error("exiting program now\n") ;
		exit(-2000) ;
    }
 
    total_allocated_RAM += n_elem * size ;
	nalloc_ptrs ++ ;
	new_cell = (malloc_cell*)malloc(sizeof(malloc_cell)) ;
	if (new_cell == NULL) {
		e_error("fatal error in memory: life functions terminated\n") ;
		exit(2001) ;
	}

	new_cell->pointer = pointer ;
	new_cell->size = n_elem * size ;
	strncpy(new_cell->filename, filename, SRCFILENAME_SZ-1) ;
	new_cell->lineno = lineno ;
	new_cell->flavour = ALLOC_CELL_RAM ;
	new_cell->next = NULL ;
	new_cell->prev = NULL ;

	/*
	 * Place new memory cell into the structure
	 */
	
	if (malloc_table == NULL) {
		malloc_table = new_cell ;
	} else {
		curr_cell = malloc_table ;
		while (curr_cell->next != NULL) {
			curr_cell = curr_cell->next ;
		}
		curr_cell->next = new_cell ;
		new_cell->prev  = curr_cell ;
	}
    return(pointer) ;
}

 
/*---------------------------------------------------------------------------
 * Function	:	ext_free_x()
 * In 		:	Pointer to de-allocate
 * Out 		:	void
 * Job		:	handles pointer free
 * Notice	:
 *--------------------------------------------------------------------------*/

void 
ext_free_x(
	void	*	pointer,
	char	*	filename,
	int			lineno
)
{
	malloc_cell	*	curr_cell ;
	int				known_flag = 0 ;

#ifdef DEBUG_ALLOC
	(void)fprintf(stderr,"<- free %p at %s (%d)\n", pointer, filename, lineno);
	(void)fflush(stderr) ;
#endif
	if (pointer == NULL) {
		e_warning("%s (%d): free requested on NULL pointer\n", filename, lineno);
		return ;
	}

	if (malloc_table == NULL) {
		e_warning("%s (%d): free requested on unallocated pointer\n",
				filename, lineno) ;
		return ;
	}

	/*
	 * Look for the requested pointer to free in our list
	 */
	curr_cell = malloc_table ;
	while (curr_cell != NULL) {
		if (curr_cell->pointer == pointer) {
			known_flag = 1 ;
			break ;
		}
		curr_cell = curr_cell->next ;
	}
	if (known_flag) {
		if (curr_cell->flavour == ALLOC_CELL_SWAP) {
			/*
			 * do a swap free: this code is actually a replicate of
			 * swap_free(), it saves a function call and a tree search
			 * to copy/paste it here.
			 */
			total_allocated_SWAP -= curr_cell->size ;
			nalloc_ptrs -- ;
			if (munmap((char*)curr_cell->pointer, curr_cell->size)!=0) {
				e_error("failed to munmap pointer %s (%d)\n", filename, lineno);
			}
			if (close(curr_cell->swapfd)==-1) {
				perror("close") ;
				e_error("closing swap file [%s]", curr_cell->swapfilename);
			}
			if (remove(curr_cell->swapfilename)!=0) {
				e_error("failed to remove file %s: errno is %d\n",
						curr_cell->swapfilename, errno) ;
				perror("remove") ;
			}
			nswapfiles -- ;
			if ((curr_cell->prev == NULL) && (curr_cell->next == NULL)) {
				malloc_table = NULL ;
			} else {
				if (curr_cell->prev == NULL) {
					/*
					 * First cell becomes malloc_table
					 */
					malloc_table = curr_cell->next ;
					malloc_table->prev = NULL ;
				} else {
					curr_cell->prev->next = curr_cell->next ;
				}
				if (curr_cell->next != NULL) {
					curr_cell->next->prev = curr_cell->prev ;
				}
			}
			free(curr_cell) ;
		} else {
			/*
			 * do a ram free: this code is actually a replicate of
			 * ram_free(), it saves a function call and a tree search
			 * to copy/paste it here.
			 */
			free(curr_cell->pointer) ;
			total_allocated_RAM -= curr_cell->size ;
			nalloc_ptrs -- ;
			if ((curr_cell->prev == NULL) && (curr_cell->next == NULL)) {
				malloc_table = NULL ;
			} else {
				if (curr_cell->prev == NULL) {
					/*
					 * First cell becomes malloc_table
					 */
					malloc_table = curr_cell->next ;
					malloc_table->prev = NULL ;
				} else {
					curr_cell->prev->next = curr_cell->next ;
				}
				if (curr_cell->next != NULL) {
					curr_cell->next->prev = curr_cell->prev ;
				}
			}
			free(curr_cell) ;
		}
	} else {
		e_warning("%s (%d) free requested on unallocated pointer\n",
				filename, lineno) ;
	}
	return ;
}


/*---------------------------------------------------------------------------
   Function	:	ram_free()
   In 		:	pointer to free, filename and lineno
   Out 		:	void
   Job		:	frees the RAM associated pointer
   Notice	:
 ---------------------------------------------------------------------------*/

void ram_free(
	void 	*	pointer,
	char	*	filename,
	int			lineno
)
{
	malloc_cell	*	curr_cell ;
	int				known_flag = 0 ;

	if (pointer == NULL) {
		e_warning("%s (%d): free requested on NULL pointer\n", filename, lineno);
		return ;
	}

	if (malloc_table == NULL) {
		e_warning("%s (%d): free requested on unallocated pointer\n",
				filename, lineno) ;
		return ;
	}

	/*
	 * Look for the requested pointer to free in our list
	 */
	curr_cell = malloc_table ;
	while (curr_cell != NULL) {
		if (curr_cell->pointer == pointer) {
			known_flag = 1 ;
			break ;
		}
		curr_cell = curr_cell->next ;
	}
	if (known_flag) {
		if (curr_cell->flavour != ALLOC_CELL_RAM) {
			e_error("trying to RAM free a non-RAM pointer: skipping free\n") ;
			return ;
		} else {
			free(curr_cell->pointer) ;
			total_allocated_RAM -= curr_cell->size ;
			nalloc_ptrs -- ;
			if ((curr_cell->prev == NULL) && (curr_cell->next == NULL)) {
				malloc_table = NULL ;
			} else {
				if (curr_cell->prev == NULL) {
					/*
					 * First cell becomes malloc_table
					 */
					malloc_table = curr_cell->next ;
					malloc_table->prev = NULL ;
				} else {
					curr_cell->prev->next = curr_cell->next ;
				}
				if (curr_cell->next != NULL) {
					curr_cell->next->prev = curr_cell->prev ;
				}
			}
			free(curr_cell) ;
		}
	} else {
		e_warning("%s (%d) free requested on unallocated pointer\n",
				filename, lineno) ;
	}
	return ;
}



/*---------------------------------------------------------------------------
   Function	:	swap_malloc_x()
   In 		:	Number of bytes to allocate in swap area
   Out 		:	pointer to allocated zone
   Job		:	Allocates memory in swap space, returns NULL if cannot
   Notice	:
   				memory should be freed with swap_free_monit()
				this function should not be called directly but through
				the swap_malloc() macro
 ---------------------------------------------------------------------------*/

void *
swap_malloc_x(
	size_t		size,
	char	*	filename,
	int			lineno
)
{
	void		*	pointer ;
	malloc_cell	*	new_cell ;
	malloc_cell	*	curr_cell ;

	char		*	fname ;
	int			 	swapfiled ;
	char		*	wbuf ;
	int				nbufs ;
	int				i ;

	if (size+total_allocated_SWAP > emem_vm_limit*ONE_MEG) {
		e_error("fatal: swap space overflow --> killing process\n") ;
		exit(-2000) ;
	}

	wbuf = (char*)calloc(emem_vm_pagesize,1) ;
	if (wbuf == NULL) {
		e_error("fatal: cannot allocate internal buffer\n") ;
		e_error("exiting program now\n") ;
		exit(-2000) ;
	}

	/* create swap file with rights: rwxrwxrwx */
	fname = get_tempfile_name() ;
 	swapfiled = open(fname, O_RDWR | O_CREAT) ;
	fchmod(swapfiled, S_IRWXU | S_IRWXG | S_IRWXO) ;

    if (size % emem_vm_pagesize == 0) {
        nbufs = size / emem_vm_pagesize ;
    } else {
        nbufs = 1 + (size / emem_vm_pagesize) ;
    }
    for (i=0 ; i<nbufs ; i++) {
        if (write(swapfiled, wbuf, emem_vm_pagesize) == -1) {
			perror("write") ;
			e_error("fatal: cannot write to swapfile: [%s] full?\n", emem_tmpdirname);
			e_error("exiting program now\n") ;
			free(wbuf) ;
			close(swapfiled) ;
			remove(fname) ;
			free(fname) ;
			exit(-2000) ;
		}
    }
	free(wbuf) ;

    /* mmap() the swap file */
    pointer = (void*)mmap(0,
                          size,
                          PROT_READ | PROT_WRITE,
                          MAP_SHARED,
                          swapfiled,
                          0) ;
    if (pointer == (void*)-1) {
        e_error("mmap failed with errno = %d\n", errno) ;
		e_error("failed to allocate SWAP %ld bytes in %s (%d)\n", 
				(long)size, filename, lineno) ;
		e_error("total SWAP allocated so far : %ld bytes\n",
		                (long)total_allocated_SWAP);
		e_error("%d swap files opened so far in %s\n", 
		                nswapfiles, emem_tmpdirname) ;
        close(swapfiled) ;
		e_error("exiting program now\n") ;
		exit(-2000) ;
    }
#ifdef DEBUG_ALLOC
	(void)fprintf(stderr,
				  "-> %p = SWAP malloc(%ld) by %s (%d) in %s/%s\n",
				  pointer, size, filename, lineno, emem_tmpdirname, fname) ;
	(void)fflush(stderr) ;
#endif

	total_allocated_SWAP += size ;	
	nalloc_ptrs ++ ;
	nswapfiles++ ;

	new_cell = (malloc_cell*)malloc(sizeof(malloc_cell)) ;
	if (new_cell == (malloc_cell*)NULL) {
		e_error("fatal error in memory: life functions terminated\n") ;
		exit(2001) ;
	}
	new_cell->pointer = pointer ;
	new_cell->size = size ;
	strncpy(new_cell->filename, filename, SRCFILENAME_SZ-1) ;
	new_cell->lineno = lineno ; 
	new_cell->flavour = ALLOC_CELL_SWAP ; 
	strncpy(new_cell->swapfilename, fname, TMPFILENAME_SZ-1) ;
	new_cell->swapfd = swapfiled ;
	new_cell->next = NULL ;
	new_cell->prev = NULL ;
		
	if (malloc_table == NULL) {
		malloc_table = new_cell ;
	} else {
		curr_cell = malloc_table ;
		while (curr_cell->next != NULL) {
			curr_cell = curr_cell->next ;
		}
		curr_cell->next = new_cell ;
		new_cell->prev  = curr_cell ;
	}
	return(pointer) ;
}

/*---------------------------------------------------------------------------
 * Function	:	swap_free_x()
 * In 		:	SWAP pointer to de-allocate
 * Out 		:	void
 * Job		:	handles swap pointer free
 * Notice	:	don't even try to understand...
 *--------------------------------------------------------------------------*/

void 
swap_free_x(
	void	*	pointer,
	char	*	filename,
	int			lineno
)
{
	malloc_cell	*	curr_cell ;
	int				known_flag = 0 ;

	if (pointer == NULL) {
		e_warning("%s (%d): free requested on NULL pointer\n", filename, lineno);
		return ;
	}

	if (malloc_table == NULL) {
		e_warning("%s (%d): free requested on unallocated pointer\n",
				filename, lineno) ;
		return ;
	}

	/*
	 * Look for the requested pointer to free in our list
	 */
	curr_cell = malloc_table ;
	while (curr_cell != NULL) {
		if (curr_cell->pointer == pointer) {
			known_flag = 1 ;
			break ;
		}
	}
	if (known_flag) {
		if (curr_cell->flavour != ALLOC_CELL_SWAP) {
			e_error("trying to RAM free a non-SWAP pointer: skipping free\n") ;
			return ;
		} else {
			total_allocated_SWAP -= curr_cell->size ;
			nalloc_ptrs -- ;
			if (munmap((char*)curr_cell->pointer, curr_cell->size)!=0) {
				e_error("failed to munmap pointer %s (%d)\n", filename, lineno);
			}
			if (close(curr_cell->swapfd)==-1) {
				perror("close") ;
				e_error("closing file [%s]", curr_cell->swapfilename);
			}
			if (close(curr_cell->swapfd)==-1) {
				perror("close") ;
				e_error("closing swap file [%s]", curr_cell->swapfilename);
			}
			if (remove(curr_cell->swapfilename)!=0) {
				e_error("failed to remove file %s: errno is %d\n",
						curr_cell->swapfilename, errno) ;
			}
			nswapfiles -- ;
			if ((curr_cell->prev == NULL) && (curr_cell->next == NULL)) {
				malloc_table = NULL ;
			} else {
				if (curr_cell->prev == NULL) {
					/*
					 * First cell becomes malloc_table
					 */
					malloc_table = curr_cell->next ;
					malloc_table->prev = NULL ;
				} else {
					curr_cell->prev->next = curr_cell->next ;
				}
				if (curr_cell->next != NULL) {
					curr_cell->next->prev = curr_cell->prev ;
				}
			}
			free(curr_cell) ;
		}
	} else {
		e_warning("%s (%d) free requested on unallocated pointer\n",
				filename, lineno) ;
	}
	return ;
}


/*---------------------------------------------------------------------------
   Function	:	ext_strdup_x()
   In 		:	char *
   Out 		:	char *
   Job		:	copy a string to a newly allocated string
   Notice	: 	returned string must be freed by free()
 ---------------------------------------------------------------------------*/

char * ext_strdup_x(char * s, char * file, int lineno)
{
	char * s2 ;

	s2 = (char*)ram_malloc(1+(int)strlen(s), file, lineno) ;
	strcpy(s2, s) ;
	return s2 ;
}



/*----------------------------------------------------------------------------
 * Function	:	get_memory_parameter()
 * In 		:	allocated character string	
 * Out 		:	size_t
 * Job		:	get info about memory handling configuration
 * Notice	:	The I/O mapping is the following:

 	Input char string				output value

	"max_ram"						maximum RAM to allocate
	"max_swap"						maximum SWAP to allocate
	"vm_page_size"					page size for virtual memory
	"total_alloc"					total allocated memory so far
	"total_alloc_ram"				total allocated RAM so far
	"total_alloc_swap"				total allocated SWAP so far

 *--------------------------------------------------------------------------*/

size_t
get_memory_parameter(char * s)
{
	if (!strcmp(s, "max_ram")) {
		return emem_ramlimit ;
	} else if (!strcmp(s, "max_swap")) {
		return emem_vm_limit ;
	} else if (!strcmp(s, "vm_page_size")) {
		return emem_vm_pagesize ;
	} else if (!strcmp(s, "total_alloc")) {
		return (size_t)(total_allocated_RAM + total_allocated_SWAP) ;
	} else if (!strcmp(s, "total_alloc_ram")) {
		return (size_t)total_allocated_RAM ;
	} else if (!strcmp(s, "total_alloc_swap")) {
		return (size_t)total_allocated_SWAP ;
	} else {
		e_error("cannot get memory parameter: [%s]\n", s) ;
		return 0 ;
	}
}

/*----------------------------------------------------------------------------
 * Function	:	set_memory_parameter()
 * In 		:	char string, size_t	
 * Out 		:	void
 * Job		:	set info about memory handling configuration
 * Notice	:	The I/O mapping is the following:

 	Input char string				output value

	"max_ram"						maximum RAM to allocate
	"max_swap"						maximum SWAP to allocate
	"vm_page_size"					page size for virtual memory

 *--------------------------------------------------------------------------*/

void
set_memory_parameter(char * s, size_t val)
{
	if (!strcmp(s, "max_ram")) {
		emem_ramlimit = val ;
	} else if (!strcmp(s, "max_swap")) {
		emem_vm_limit = val ;
	} else if (!strcmp(s, "vm_page_size")) {
		emem_vm_pagesize = val ;
	} else {
		e_error("cannot set memory parameter: [%s]\n", s) ;
	}
	return ;
}


/*---------------------------------------------------------------------------
   Function	:	print_memory_parameters()
   In 		:	void
   Out 		:	void, messages to stderr
   Job		:	print out current memory configuration
   Notice	:
 ---------------------------------------------------------------------------*/

void print_memory_parameters(void)
{
	fprintf(stderr, "MEM: ramlimit  %ld\n", (long)emem_ramlimit) ;
	fprintf(stderr, "MEM: vmlimit   %ld\n", (long)emem_vm_limit) ;
	fprintf(stderr, "MEM: pagesize  %ld\n", (long)emem_vm_pagesize) ;
	fprintf(stderr, "MEM: tmpdir    [%s]\n", emem_tmpdirname) ;
}



/*---------------------------------------------------------------------------
 * Function	:	memory_status()
 * In 		:	void
 * Out 		:	prints out memory status on stderr
 * Notice	:
 *--------------------------------------------------------------------------*/
void	
memory_status(void)
{
	malloc_cell	*curr_cell ;

	if ((total_allocated_RAM!=0) ||
		(total_allocated_SWAP!=0) ||
		(nswapfiles>0) || 
		(nalloc_ptrs>0)) {
		(void)fprintf(stderr, "**** MEMORY STATUS\n") ;
		(void)fprintf(stderr, "memory in use  : %ld bytes\n",
				(long)total_allocated_RAM + total_allocated_SWAP) ;
		(void)fprintf(stderr, "RAM            : %d bytes\n",
				total_allocated_RAM);
		(void)fprintf(stderr, "SWAP           : %ld bytes\n",
				(long)total_allocated_SWAP);
		(void)fprintf(stderr, "open swapfiles : %d\n", nswapfiles) ;
		
		curr_cell = malloc_table ;
		while (curr_cell != NULL) {

			if (curr_cell->flavour == ALLOC_CELL_RAM) {
				(void)fprintf(stderr, "[RAM] ") ;
			} else {
				(void)fprintf(stderr, "[SWP] ") ;
			}

			(void)fprintf(stderr, "po=%p size=%ld allocated in %s (%d) ",
					curr_cell->pointer,
					(long)curr_cell->size,
					curr_cell->filename,
					(int) curr_cell->lineno) ;

			if (curr_cell->flavour == ALLOC_CELL_SWAP) {
				(void)fprintf(stderr, "in swapfile %s",
						curr_cell->swapfilename) ;
			}
			(void)fprintf(stderr, "\n") ;
			(void)fflush(stderr) ;

			curr_cell = curr_cell->next ;
		}
	}
	return ;
}


/*----------------------------------------------------------------------------
 * Function	:	cleanup_mess()
 * In 		:	void	
 * Out 		:	void	
 * Job		:	clean up temporary files created by the current process
 * Notice	:	
 *--------------------------------------------------------------------------*/

void cleanup_mess(void)
{
	pid_t	pid ;
	char	name[512] ;
	int		i ;

	pid = getpid() ;
	if (pid == (pid_t)-1) {
		e_error("cannot get current process ID: no cleanup of tmp files\n") ;
		return ;
	}
	
	for (i=1 ; i<=file_reg ; i++) {
		(void)sprintf(name, "%s/vmswap_%05ld_%05x",
					  emem_tmpdirname,
					  (long)pid,
					  (unsigned int)i) ;
		if (remove(name) == 0) {
			fprintf(stderr, "\n*** cleaned up tmp file [%s]", name) ;
		}
	}
	fprintf(stderr, "\n\n") ;
	fflush(stderr) ;
}

void cleanup_mess_wrapper(void)
{ 
	fprintf(stderr, "\n\n**** user interrupted process %ld ****\n",
			(long)getpid()) ;
	exit(1) ;
}


/*----------------------------------------------------------------------------
 * Function	:	signal_catch()
 * In 		:	signal to catch, function to call in this case
 * Out 		:	int 0 if worked, -1 otherwise
 * Job		:	catch an interrupt and launch the appropriate function
 * Notice	:	directly taken from MIDAS
 *--------------------------------------------------------------------------*/
typedef void (*func)(int);
int signal_catch(int sig, void (*f)())
{
	struct sigaction	act,
						oact;
	
	/*
	 * Simply install the given function f() as signal handler for
	 * interrupts
	 */
	act.sa_handler = (func)f;
	sigemptyset(&act.sa_mask);
	act.sa_flags = 0;
	if (sigaction(sig,&act,&oact) != 0) return -1;
	return 0 ;
}


/*---------------------------------------------------------------------------
   Function	:	fatal_signal_handler()
   In 		:	void
   Out 		:	void
   Job		:	prints out a message to stderr and call cleanup
   Notice	:	called from SIGSEGV, SIGBUS, SIGXCPU, SIGXFSZ
 ---------------------------------------------------------------------------*/

void fatal_signal_handler(void)
{
	fprintf(stderr, "\n\n") ;
	fprintf(stderr, "**** fatal error in pid %ld\n", (long)getpid()) ;
	fprintf(stderr, "**** Segmentation fault or Bus Error\n") ;
	exit(1) ;
}


 



/*---------------------------------------------------------------------------
   Function	:	mmap_open()
   In 		:	filename
   Out 		:	pointer to allocated mmap_file structure
   Job		:	mmap the totality of a given file
   Notice	:	should only be closed by mmap_close()
 ---------------------------------------------------------------------------*/
size_t mem_filesize(char *filename);
mmap_file *
mmap_open(char * filename)
{
	mmap_file	*	mm ;
	char		*	buf ;
	long			size ;
	int				fd ;

	if ((size = mem_filesize(filename)) < 1) return NULL ;
	if ((fd = open(filename, O_RDONLY)) == -1) {
		e_error("cannot open file %s: aborting mmapping", filename) ;
		return NULL ;
	}
	buf = (char*)mmap(0, size, PROT_READ, MAP_SHARED, fd, 0) ;
	if (buf==(char*)-1) { perror("mmap") ; close(fd) ; return NULL ; }
	close(fd) ;

	mm = (mmap_file*)malloc(sizeof(mmap_file)) ;
	mm->fd 	 = fd ;
	mm->buf  = buf ;
	mm->size = size ;
	(void)strcpy(mm->filename, filename);
	return mm ;
}


/*---------------------------------------------------------------------------
   Function	:	mmap_close()
   In 		:	allocated mmapped_file struct pointer
   Out 		:	void
   Job		:	munmap and close the file
   Notice	:
 ---------------------------------------------------------------------------*/

void mmap_close(mmap_file * mm)
{
	if (mm==NULL) return ;
	if (munmap((char*)mm->buf, mm->size)!=0) perror("munmap") ;
	free(mm);
	return ;
}



/*---------------------------------------------------------------------------
   Function	:	dump_stack()
   In 		:	void
   Out 		:	void, messages on stderr
   Job		:	dump the current stack to stderr, using dbx
   Notice	:	Probably Solaris specific, with little efforts portable
   				other Unixes. From the comp.unix.programmer FAQ.
 ---------------------------------------------------------------------------*/

void dump_stack(void)
{
	char s[160];
	sprintf(s, "/bin/echo 'where\ndetach' | dbx - %d", (int) getpid());
	system(s);
	return;
}

#ifdef EXPORT_MEMORY_H 
size_t mem_filesize(char *filename)
{
    size_t size ;
    struct stat fileinfo ;
    /* POSIX compliant  */
    if (stat(filename, &fileinfo) != 0) {
        size = (size_t)0 ;
    } else {
        size = (size_t)fileinfo.st_size ;
    }
    return size ;
}
#endif
}
/*--------------------------------------------------------------------------*/
#endif



