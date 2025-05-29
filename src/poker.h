#pragma once

#include "common.h"

namespace Poker {
    char const SUIT_SYMBOLS[] = {'C', 'S', 'D', 'H'};
    char const RANK_SYMBOLS[] = {'2', '3', '4', '5', '6', '7', '8', '9', 'T', 'J', 'Q', 'K', 'A'};

    struct Card {
        int rank = 0;
        int suit = 0;

        Card() {}

        template<typename Cont>
        Card(Cont const& str) {
            while(RANK_SYMBOLS[rank] != str[0]) checkval(++rank, assert_(_ < 13));
            while(SUIT_SYMBOLS[suit] != str[1]) checkval(++suit, assert_(_ < 4));
        }

        bool operator==(Card const& other) const& {
            return (rank == other.rank) && (suit == other.suit);
        }
    };

    template<typename OStrm>
    OStrm& operator<<(OStrm& strm, Card const& card) {
        return strm << RANK_SYMBOLS[card.rank] << SUIT_SYMBOLS[card.suit];
    }

    using Hand = std::array<Card, 5>;

    enum class HandRank {
        HIGH_CARD, PAIR, TWO_PAIRS, TRIPLE, STRAIGHT, FLUSH, FULL_HOUSE, QUADRUPLE, STRAIGHT_FLUSH
    };

    auto Strength(Hand const& h) {
        assert_([&](){
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    if (i == j) continue;
                    if (h[i] == h[j]) return false;
                }
            }
            return true;
        }())
        std::array<int, 4> suits;
        std::fill(suits.begin(), suits.end(), 0);
        std::array<int, 13> ranks;
        std::fill(ranks.begin(), ranks.end(), 0);
        for (auto const c : h) {
            ++suits[c.suit]; ++ranks[c.rank];
        }
        bool flush = std::find(suits.begin(), suits.end(), 5) != suits.end();
        int straight_strength = -1;
        for (int i = -1; i < 13; ++i) {
            if (i - straight_strength == 5) {
                break;
            } else if (ranks[i == -1 ? 12 : i] != 1) {
                straight_strength = i + 1;
            }
        }
        if (flush && straight_strength <= 8/*T*/) {
            return std::make_pair(HandRank::STRAIGHT_FLUSH, straight_strength);
        }
        std::array<std::vector<int>, 5> mult;
        for (int r = 12; r >= 0; --r) {
            mult[ranks[r]].push_back(r);
        }

        int total_strength = 0;
        for (int m = 4; m > 0; --m) {
            for (auto const r : mult[m]) {
                total_strength *= 13;
                total_strength += r;
            }
        }

        if (mult[4].size() == 1) {
            return std::make_pair(HandRank::QUADRUPLE, total_strength);
        }
        if (mult[3].size() == 1 && mult[2].size() == 1) {
            return std::make_pair(HandRank::FULL_HOUSE, total_strength);
        }
        if (flush) {
            return std::make_pair(HandRank::FLUSH, total_strength);
        }
        if (straight_strength <= 8/*T*/) {
            return std::make_pair(HandRank::STRAIGHT, straight_strength);
        }
        if (mult[3].size() == 1) {
            return std::make_pair(HandRank::TRIPLE, total_strength);
        }
        if (mult[2].size() == 2) {
            return std::make_pair(HandRank::TWO_PAIRS, total_strength);
        }
        if (mult[2].size() == 1) {
            return std::make_pair(HandRank::PAIR, total_strength);
        }
        return std::make_pair(HandRank::HIGH_CARD, total_strength);
    }
}
