#pragma once

#include "common.h"

template<typename Perm>
bool AdvancePermutation(Perm& perm) {
    for (int i = perm.size() - 1; i > 0; --i) {
        if (perm[i] > perm[i - 1]) {
            int j;
            for (j = i + 1; j < perm.size(); ++j) {
                if (perm[j] <= perm[i - 1]) break;
            }
            std::swap(perm[j - 1], perm[i - 1]);
            std::reverse(perm.begin() + i, perm.end());
            return true;
        }
    }
    return false;
}

struct Blank {};

template<typename TMin = Blank, typename TMax = Blank, typename TLen = Blank>
void WithAllPartitions(auto const& fn, auto const& sum, TMin const& min_el = {}, TMax const& max_el = {}, TLen const& max_len = {}) {
    
}