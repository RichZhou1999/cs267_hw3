#pragma once
#include <math.h>
#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {
    // serial implementation
	/*
	std::vector<kmer_pair> data;
    std::vector<int> used;
	*/

	// arrays of global ptrs that point to arrays allocated in the shared
	// segment on each rank
	std::vector< upcxx::global_ptr<kmer_pair> > data;
	std::vector< upcxx::global_ptr<int> >       used;

    size_t my_size;
	size_t arr_size;
	size_t vsize;
    size_t size() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Helper functions

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);

	// helper functions
	uint64_t get_rank_from_slot(uint64_t slot);
	uint64_t get_ind_from_slot(uint64_t slot);
};

HashMap::HashMap(size_t size) {
    my_size = size;
	arr_size = (size_t)ceil(size / upcxx::rank_n() );
   	vsize = size_t(upcxx::rank_n());
	data.resize(vsize);
    used.resize(vsize, 0);
	// init arrays in the vectors
	for (int i = 0; i < vsize; i++){
		data[i] = upcxx::new_array<kmer_pair>(arr_size);
		used[i] = upcxx::new_array<int>(arr_size);
	}
}

bool HashMap::insert(const kmer_pair& kmer) {
    int rank  = upcxx::rank_me();
	int nrank = upcxx::rank_n();
	uint64_t hash = kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
        // at slot == hash do we have a value? slot cannot exceed size
		// if yes: increment probe
		uint64_t slot = (hash + probe++) % size();
		// succes true if slot is empty
        success = request_slot(slot);
        if (success) {
			// if avail slot found: write kmer to data at slot
            write_slot(slot, kmer);
        }
    } while (!success && probe < size());
    return success;
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
		// at slot == hash do we have a value?
        // if yes: slot_used returns true
        uint64_t slot = (hash + probe++) % size();
        if (slot_used(slot)) {
            val_kmer = read_slot(slot);
			// if value found success is true
            if (val_kmer.kmer == key_kmer) {
                success = true;
            }
        }
    } while (!success && probe < size());
    return success;
}

uint64_t HashMap::get_rank_from_slot(uint64_t slot){
	uint64_t rank;
	rank = uint64_t( slot/arr_size );
	return rank;
}

uint64_t HashMap::get_ind_from_slot(uint64_t slot){
    uint64_t rank;
    rank = uint64_t( slot/arr_size );
    uint64_t ind = uint64_t( slot - rank*arr_size );
	return ind;
}

bool HashMap::slot_used(uint64_t slot) { 
	return used[slot] != 0; 
}

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { 
	upcxx::global_ptr<kmer_pair> kmer_glptr = data[slot];
	// downcasting
	UPCXX_ASSERT( kmer_glptr.is_local() );
	kmer_pair* kmer_loptr = kmer_glptr.local();
	*kmer_loptr = kmer;
}

kmer_pair HashMap::read_slot(uint64_t slot) { 
	upcxx::global_ptr<kmer_pair> kmer_glptr = data[slot];
    // downcasting
    UPCXX_ASSERT( kmer_glptr.is_local() );
    kmer_pair* kmer_loptr = kmer_glptr.local();
	return *kmer_loptr; 
}

bool HashMap::request_slot(uint64_t slot) {
    upcxx::global_ptr<int> used_slot_glptr = used[slot];
    // downcasting
    UPCXX_ASSERT( used_slot_glptr.is_local() );
    int used_slot = *used_slot_glptr.local();
	if ( used_slot != 0) {
        return false;
    } else {
        used_slot = 1;
        return true;
    }
}

size_t HashMap::size() const noexcept { return my_size; }

