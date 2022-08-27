//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// lru_replacer.cpp
//
// Identification: src/buffer/lru_replacer.cpp
//
// Copyright (c) 2015-2019, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "buffer/lru_replacer.h"

namespace bustub {

LRUReplacer::LRUReplacer(size_t num_pages) 
{
    capcity = num_pages;
    size = 0;
    lists.clear();
    hash.clear();
}

LRUReplacer::~LRUReplacer() = default;

auto LRUReplacer::Victim(frame_id_t *frame_id) -> bool 
{
    if (lists.empty()) return false;
    mtx.lock();
    frame_id_t last_frame = lists.back();
    lists.pop_back();
    hash.erase(last_frame);
    mtx.unlock();
    *frame_id = last_frame;
    return true;
}

void LRUReplacer::Pin(frame_id_t frame_id) 
{
    if (!hash.count(frame_id)) return;
    mtx.lock();
    lists.erase(hash[frame_id]);
    hash.erase(frame_id);
    mtx.unlock();
}

void LRUReplacer::Unpin(frame_id_t frame_id) 
{
    if (hash.count(frame_id)) return;
    mtx.lock();
    if (lists.size() == capcity)
        lists.pop_back();
    lists.push_front(frame_id);
    std::list<frame_id_t>::iterator it = lists.begin();
    hash.insert({frame_id, it});
    mtx.unlock();
}

auto LRUReplacer::Size() -> size_t { return size = lists.size();}

}  // namespace bustub
