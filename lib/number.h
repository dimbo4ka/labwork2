#pragma once
#include <cinttypes>
#include <iostream>


struct uint239_t {
    uint8_t data[35]{};
};

static_assert(sizeof(uint239_t) == 35, "Size of uint239_t must be no higher than 35 bytes");

uint239_t FromInt(uint32_t value, uint64_t shift);

uint239_t FromString(const char* str, uint64_t shift);

uint239_t operator+(const uint239_t& lhs, const uint239_t& rhs);

uint239_t operator-(const uint239_t& lhs, const uint239_t& rhs);

uint239_t operator*(const uint239_t& lhs, const uint239_t& rhs);

uint239_t operator/(const uint239_t& lhs, const uint239_t& rhs);

bool operator==(const uint239_t& lhs, const uint239_t& rhs);

bool operator!=(const uint239_t& lhs, const uint239_t& rhs);

bool operator<(const uint239_t& lhs, const uint239_t& rhs);

bool operator>(const uint239_t& lhs, const uint239_t& rhs);

bool operator<=(const uint239_t& lhs, const uint239_t& rhs);

bool operator>=(const uint239_t& lhs, const uint239_t& rhs);

std::ostream& operator<<(std::ostream& stream, const uint239_t& value);

uint64_t GetShift(const uint239_t& value);