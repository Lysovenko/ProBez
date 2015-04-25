/*
 *      (C) Serhii Lysovenko <lisovenko.s at the Gmail>
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 3 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */

#include <stdlib.h>
#include <string.h>
#include "interf.h"
static void **TotalDataForRequests;

void
init_requests ()
{
  TotalDataForRequests = malloc (sizeof (void *) * N_requests);
  memset (TotalDataForRequests, 0, sizeof (void *) * N_requests);
}

void *
get_request (int n)
{
  return TotalDataForRequests[n];
}

void
set_request (int n, void *what)
{
  TotalDataForRequests[n] = what;
}

void
unset_request (int n)
{
  TotalDataForRequests[n] = NULL;
}
