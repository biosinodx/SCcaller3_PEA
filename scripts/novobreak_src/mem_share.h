/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef __MEM_SHARE_RJ_H
#define __MEM_SHARE_RJ_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#ifndef OBJ_DESC_MAX_CHILD
#define OBJ_DESC_MAX_CHILD 64
#endif

static inline size_t mem_size_round(size_t size){ return (size + 7) & 0xFFFFFFFFFFFFFFF8LLU; }

static inline uint8_t mem_size_gap(size_t size){ return (size & 0x07U)? 8 - (size & 0x07U) : 0; }

static inline size_t mem_dump(void *mem, size_t len, FILE *out){
	size_t size;
	uint8_t i, v;
	size = mem_size_round(len);
	if(out){
		fwrite(mem, 1, len, out);
		v = 0;
		for(i=0;i<mem_size_gap(len);i++) fwrite(&v, 1, 1, out);
	}
	return size;
}

typedef size_t (*mem_array_count)(void *obj, int idx);

struct obj_desc_t;

typedef struct obj_desc_t {
	size_t size;
	int n_child; // <= OBJ_DESC_MAX_CHILD
	off_t addr[OBJ_DESC_MAX_CHILD];
	const struct obj_desc_t *desc[OBJ_DESC_MAX_CHILD];
	mem_array_count cnt;
} obj_desc_t;

static inline size_t mem_size_obj(void *obj, const obj_desc_t *desc, size_t size, size_t cnt){
	size_t m;
	void *ref;
	int i;
	if(desc == NULL) return size;
	if(obj == NULL) return size;
	size += mem_size_round(cnt * desc->size);
	if(desc->n_child == 0) return size;
	for(m=0;m<cnt;m++){
		ref = obj + m * desc->size;
		for(i=0;i<desc->n_child;i++){
			size += mem_size_obj(*((void**)(ref + desc->addr[i])), desc->desc[i], 0, desc->cnt? desc->cnt(ref, i) : 1);
		}
	}
	return size;
}

static inline size_t mem_dump_obj(void *obj, const obj_desc_t *desc, size_t offset, size_t cnt, FILE *out, int mem_free){
	void *ref;
	size_t size, m;
	int i;
	if(obj == NULL) return offset;
	size = mem_dump(obj, cnt * desc->size, out) + offset;
	if(desc->n_child == 0) return size;
	for(m=0;m<cnt;m++){
		ref = obj + desc->size * m;
		for(i=0;i<desc->n_child;i++){
			size = mem_dump_obj(*((void**)(ref + desc->addr[i])), desc->desc[i], size, desc->cnt? desc->cnt(ref, i) : 1, out, mem_free);
		}
	}
	if(mem_free) free(obj);
	return size;
}

static inline size_t mem_dump_obj_file(void *obj, const obj_desc_t *desc, size_t cnt, size_t aux_data, FILE *out){
	size_t size;
	if(desc == NULL) return 0;
	size = mem_size_obj(obj, desc, 0, cnt);
	fwrite(&size, sizeof(size_t), 1, out);
	fwrite(&cnt, sizeof(size_t), 1, out);
	fwrite(&aux_data, sizeof(size_t), 1, out);
	size = 3 * sizeof(size_t);
	size = mem_dump_obj(obj, desc, size, cnt, out, 0);
	fflush(out);
	return size;
}

static inline size_t mem_dump_free_obj_file(void *obj, const obj_desc_t *desc, size_t cnt, size_t aux_data, FILE *out){
	size_t size;
	if(desc == NULL) return 0;
	size = mem_size_obj(obj, desc, 0, cnt);
	if(out){
		fwrite(&size, sizeof(size_t), 1, out);
		fwrite(&cnt, sizeof(size_t), 1, out);
		fwrite(&aux_data, sizeof(size_t), 1, out);
	}
	size = 3 * sizeof(size_t);
	size = mem_dump_obj(obj, desc, size, cnt, out, 1);
	fflush(out);
	return size;
}

static inline size_t mem_load_obj(void *obj, const obj_desc_t *desc, size_t cnt){
	size_t size, m;
	int i;
	void *ref, **ptr;
	if(desc == NULL) return 0;
	size = mem_size_round(desc->size * cnt);
	if(desc->n_child == 0) return size;
	for(m=0;m<cnt;m++){
		ref = obj + desc->size * m;
		for(i=0;i<desc->n_child;i++){
			ptr = (void**)(ref + desc->addr[i]);
			if(*ptr == NULL) continue;
			*ptr = obj + size;
			size += mem_load_obj(*ptr, desc->desc[i], desc->cnt? desc->cnt(ref, i) : 1);
		}
	}
	return size;
}

static char *mem_share_locks = NULL;
static int mem_share_lock_size = 0;

static inline void cleanup_mem_share_file_locks(){
	int off;
	off = 0;
	while(off < mem_share_lock_size){
		unlink(mem_share_locks + off);
		off += strlen(mem_share_locks + off) + 1;
	}
	if(mem_share_locks) free(mem_share_locks);
	mem_share_lock_size = 0;
	mem_share_locks = NULL;
}

#ifndef sighandler_t
typedef void (*sighandler_t)(int sig);
#endif
static sighandler_t sig_term = SIG_IGN;
static sighandler_t sig_int  = SIG_IGN;
static sighandler_t sig_hup  = SIG_IGN;
static sighandler_t sig_kill = SIG_IGN;
static volatile sig_atomic_t cleanup_mem_share_in_progress = 0;

static inline void sig_cleanup_mem_share_file_locks(int sig){
	if(cleanup_mem_share_in_progress) raise(sig);
	cleanup_mem_share_in_progress = 1;
	cleanup_mem_share_file_locks();
	signal(SIGTERM, sig_term);
	signal(SIGINT , sig_int);
	signal(SIGHUP, sig_hup);
	signal(SIGKILL, sig_kill);
	raise(sig);
}

static inline void register_mem_share_file_lock(char *file){
	int len;
	if(mem_share_lock_size == 0){
		if((sig_term = signal(SIGTERM, sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGTERM, SIG_IGN);
		if((sig_int  = signal(SIGINT , sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGINT , SIG_IGN);
		if((sig_hup  = signal(SIGHUP , sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGHUP , SIG_IGN);
		if((sig_kill = signal(SIGKILL, sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGKILL, SIG_IGN);
		atexit(cleanup_mem_share_file_locks);
	}
	len = strlen(file);
	mem_share_locks = realloc(mem_share_locks, mem_share_lock_size + len + 1);
	strcpy(mem_share_locks + mem_share_lock_size, file);
	mem_share_lock_size += len + 1;
}

static inline void* mem_load_obj_file_core(const obj_desc_t *desc, char *path, size_t *cnt, size_t *aux_data){
	void *mem;
	size_t size, psize, i, tot;
	char *lock;
	FILE *file;
	int fd;
	if(desc == NULL) return NULL;
	if((file = fopen(path, "r+")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
	}
	fread(&size, sizeof(size_t), 1, file);
	fread(cnt, sizeof(size_t), 1, file);
	fread(aux_data, sizeof(size_t), 1, file);
	fd = fileno(file);
	psize = getpagesize();
	mem = mmap(0, (size + 3 * sizeof(size_t) + psize - 1) / psize * psize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if(mem == MAP_FAILED){
		perror("Cannot mmap");
		abort();
	}
	fclose(file);
	mem_load_obj(mem + 3 * sizeof(size_t), desc, *cnt);
	lock = alloca(strlen(path) + 32);
	sprintf(lock, "%s.mem_share.%ld", path, gethostid());
	if((file = fopen(lock, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", lock, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
	}
	fwrite(&mem, sizeof(void*), 1, file);
	fclose(file);
	tot = 0;
	for(i=0;i<size;i++){
		tot += ((char*)(mem + 3 * sizeof(size_t)))[i];
	}
	if(tot == 0){
		fprintf(stderr, "-- ruanjue --\n"); fflush(stderr);
	}
	fprintf(stderr, "-- mem_share '%s' at %p --\n", path, mem); fflush(stderr);
	register_mem_share_file_lock(lock);
	return mem + 3 * sizeof(size_t);
}

static inline void* mem_load_obj_file(const obj_desc_t *desc, char *path, size_t *cnt, size_t *aux_data){
	char *lock;
	void *addr, *mem;
	size_t size, psize;
	FILE *file;
	int fd;
	lock = alloca(strlen(path) + 32);
	sprintf(lock, "%s.mem_share.%ld", path, gethostid());
	if((file = fopen(lock, "r")) == NULL){
		return mem_load_obj_file_core(desc, path, cnt, aux_data);
	}
	fread(&addr, sizeof(void*), 1, file);
	fclose(file);
	if((file = fopen(path, "r")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
	}
	fd = fileno(file);
	fread(&size, sizeof(void*), 1, file);
	fread(cnt, sizeof(void*), 1, file);
	fread(aux_data, sizeof(void*), 1, file);
	psize = getpagesize();
	mem = mmap(addr, (size + 3 * sizeof(size_t) + psize - 1) / psize * psize, PROT_READ, MAP_SHARED | MAP_FIXED, fd, 0);
	fclose(file);
	if(mem == MAP_FAILED){
		perror("Cannot map shared object");
		return mem_load_obj_file_core(desc, path, cnt, aux_data);
	}
	fprintf(stderr, "-- mem_map '%s' at %p --\n", path, addr); fflush(stderr);
	return mem + 3 * sizeof(size_t);
}

static const struct obj_desc_t OBJ_DESC_ARRAY = {1, 0, {}, {}, NULL};

/**
 * To define a obj_desc_t for a array of pointers, such as u32hash*
 *  u32hash **hashs];
 *  static struct obj_desc_t u32hashv_obj_desc = {sizeof(u32hash*), 1, {0}, {(obj_desc_t*)&u32hash_obj_desc}, NULL};
 * */

#endif
