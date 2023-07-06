/*
 * header fiel from Advanced Programming in the UNIX Environment 
 */
#ifndef _APUE_H
#define _APUE_H

#define _POSIX_C_SOURCE 200809L

#if defined(SOLARIS)
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 700
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/termios.h>

#if defined(MACOS) || !defined(TIOCGWINSZ)
#include <sys/ioctl.h>
#endif

#include <stdio.h>  /* for convenience */
#include <stdlib.h> /* for convenience */
#include <stddef.h> /* for offsetof */
#include <string.h> /* for convenience */
#include <unistd.h> /* for convenience */
#include <signal.h> /* for  SIG_ERR */

#define MAXLINE 4096  /* MAX line length */

/*
 * Default file access permissions for new files.
 */
#define FILE_MODE (S_IRUSE | S_IWUSR | S_RGRP | S_IROTH)

/*
 * Defualt permmsion for new directiories.
 */
#define DIR_MODE (FILE_MODE | S_IXUSR | S_IXGRP | S_IXOTH) 



#endif
