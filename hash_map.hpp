#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {
    upcxx::atomic_domain<int> ad = upcxx::atomic_domain<int>({upcxx::atomic_op::compare_exchange,upcxx::atomic_op::load});

    std::vector <upcxx::global_ptr<kmer_pair>> global_data;
    std::vector <upcxx::global_ptr<int>> global_used;

    size_t global_size;
    size_t my_size;
    int my_rank;
    int number_of_ranks;
    int start_position;
    int offset;

    size_t size() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    int compare_exchange(uint64_t);
    // bool slot_used(uint64_t slot);
    bool slot_used(uint64_t pos);//, upcxx::global_ptr<int> tmp);

    // Helper functions

    // Request a slot or check if it's already used.
    // bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t size) {
    this->global_size = size;
    number_of_ranks = upcxx::rank_n();
    my_rank =  upcxx::rank_me();
    offset = (global_size + number_of_ranks - 1) / number_of_ranks;
    start_position = offset * my_rank;
    my_size = offset;
    if (my_rank == number_of_ranks - 1) {
        my_size = global_size - (number_of_ranks - 1) * offset;
    }
    global_data = std::vector<upcxx::global_ptr<kmer_pair>> (number_of_ranks);
    global_used = std::vector<upcxx::global_ptr<int>> (number_of_ranks);

    upcxx::global_ptr<kmer_pair> local_data = upcxx::new_array<kmer_pair>(my_size);
    upcxx::global_ptr<int> local_used = upcxx::new_array<int>(my_size);
    upcxx::future<upcxx::global_ptr<kmer_pair>> * data_fut = new upcxx::future<upcxx::global_ptr<kmer_pair>>[number_of_ranks];
    upcxx::future<upcxx::global_ptr<int>> * used_fut = new upcxx::future<upcxx::global_ptr<int>>[number_of_ranks];
    for (int i = 0; i < upcxx::rank_n(); i++) {
        data_fut[i] = upcxx::broadcast(local_data, i);
        used_fut[i] = upcxx::broadcast(local_used, i);
    }
    for (int i = 0; i < upcxx::rank_n(); i++) {
        global_data[i] = data_fut[i].wait();
        global_used[i] = used_fut[i].wait();
    }
}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    int pos = hash % global_size;
    int start_pos = pos;
    while(compare_exchange(pos)!=0){
        pos += 1;
        pos = pos % global_size;
        if(start_pos == pos){
            return false;
        }
    }
    upcxx::rput(kmer, global_data[pos/offset]+(pos%offset)).wait();
    return true;

}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    int pos = hash % global_size;
    int start_pos = pos;

    do{
        if(slot_used(pos)){
            return false;
        }
        kmer_pair retreived_kmer = upcxx::rget(global_data[pos/offset]+(pos%offset)).wait();
        if (retreived_kmer.kmer == key_kmer){
            val_kmer = retreived_kmer;
            return true;
        }
        pos += 1;
        pos = pos % global_size;
    }while (pos != start_pos);
    return false;

}

bool HashMap::slot_used(uint64_t pos) {
    return ad.load(global_used[pos/offset] + (pos%offset), std::memory_order_relaxed).wait() == 0; }


int HashMap::compare_exchange(uint64_t pos) {
    return ad.compare_exchange(global_used[pos/offset] + (pos%offset),
                               0, 1, std::memory_order_relaxed).wait();
}



size_t HashMap::size() const noexcept { return my_size; }



