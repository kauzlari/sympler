/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2013, 
 * David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
 * and others authors stated in the AUTHORS file in the top-level 
 * source directory.
 *
 * SYMPLER is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SYMPLER is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SYMPLER.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Please cite the research papers on SYMPLER in your own publications. 
 * Check out the PUBLICATIONS file in the top-level source directory.
 *
 * You are very welcome to contribute extensions to the code. Please do 
 * so by making a pull request on https://github.com/kauzlari/sympler
 * 
 */


#include "threads.h"

#include "general.h"

#ifdef ENABLE_PTHREADS

pthread_mutexattr_t g_mutex_attr;


void *RunThread(thread_info_t *t)                              
{                                                         
  ((BasicThread*) t->handler)->operator()(t->no);     
                                                          
  pthread_exit(NULL);                                     
}


void ForkThreads(BasicThread *handler)
{
  ForkThreadsWithFunction((thread_func_t) RunThread, handler);
}


void ForkThreadsWithFunction(thread_func_t function, void *handler)
{
  thread_info_t *thread = new thread_info_t[global::n_threads];

  /*
  MSG_DEBUG
    ("--- ForkThreadsWithFunction",
     "Forking " << global::n_threads << " threads.");
  */

  for (int t = 0; t < global::n_threads; t++) {
    thread[t].no = t;
    thread[t].handler = handler;
    pthread_create(&thread[t].id, NULL, function, &thread[t]);
  }

  for (int t = 0; t < global::n_threads; t++) {
    pthread_join(thread[t].id, NULL);
  }

  /*
  MSG_DEBUG
    ("--- ForkThreadsWithFunction",
     "Synchronized.");
  */
  delete thread;
}

#endif
