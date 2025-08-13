#include "number.h"

const uint8_t kNonServiceBitsInByte = 7;
const uint8_t kPaddingBitsOfType = 6;
const uint8_t kSignificantBitsOfType = 239;
const uint8_t kUint239Size = 35;
const uint64_t kShiftMaxValue = (1ll << 35) - 1ll;
const uint239_t kUint239Zero = FromInt(0, 0); 

uint8_t GetBitValue(const uint239_t& number, uint8_t byte, uint8_t bit) {
    return (number.data[byte] >> bit) & 1;
}

uint8_t GetStartBit(uint64_t shift) {
    return (
            (kSignificantBitsOfType + kPaddingBitsOfType + shift)
             % (kSignificantBitsOfType + kPaddingBitsOfType)
           )
            % kNonServiceBitsInByte;
}

uint8_t GetStartByte(uint64_t shift) {
    return kUint239Size - 1 
           - (
              (kSignificantBitsOfType + kPaddingBitsOfType + shift) 
               % (kSignificantBitsOfType + kPaddingBitsOfType)
             )
              / kNonServiceBitsInByte;
}

void MoveToNextBit(uint8_t& current_bit, uint8_t& current_byte) {
    current_bit = (current_bit + 1) % kNonServiceBitsInByte;
    if (current_bit != 0)
        return;

    current_byte = (kUint239Size + current_byte - 1) % kUint239Size;
}

void MoveNumberToLeft(uint239_t& number, uint64_t moves) {
    uint239_t result = FromInt(0, GetShift(number));
    
    uint8_t current_number_bit = GetStartBit(GetShift(number));
    uint8_t current_number_byte = GetStartByte(GetShift(number));
    uint8_t current_result_bit = GetStartBit(GetShift(number) + moves);
    uint8_t current_result_byte = GetStartByte(GetShift(number) + moves);

    for (size_t i = 0; i < kSignificantBitsOfType - moves; ++i) {
        uint8_t number_bit = GetBitValue(
                                         number
                                         , current_number_byte
                                         , current_number_bit
                                        );
        result.data[current_result_byte] += (number_bit << current_result_bit);

        MoveToNextBit(current_number_bit, current_number_byte);
        MoveToNextBit(current_result_bit, current_result_byte);
    }

    number = result;
}

void ChangeShift(uint239_t& number, uint64_t shift) {
    uint239_t value = FromInt(0, shift);

    uint8_t current_number_bit = GetStartBit(GetShift(number));
    uint8_t current_number_byte = GetStartByte(GetShift(number));
    uint8_t current_value_bit = GetStartBit(shift);
    uint8_t current_value_byte = GetStartByte(shift);

    for(size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t number_bit = GetBitValue(
                                         number
                                         , current_number_byte
                                         , current_number_bit
                                        );
        value.data[current_value_byte] += (number_bit) << current_value_bit;

        MoveToNextBit(current_number_bit, current_number_byte);
        MoveToNextBit(current_value_bit, current_value_byte);
    }
    
    number = value;
}

void SetShift(uint239_t& number, uint64_t shift) {
    for (size_t current_byte = 0; current_byte < kUint239Size; ++current_byte) {
        number.data[kUint239Size - 1 - current_byte]
            += ((shift >> current_byte) & 1) << 7;
    }
}

bool operator<(const uint239_t& lhs, const uint239_t& rhs) {
    for (size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t current_lhs_bit
            = GetStartBit(
                          GetShift(lhs) - i - 1 - kPaddingBitsOfType
                         );
        uint8_t current_lhs_byte
            = GetStartByte(
                           GetShift(lhs) - i - 1 - kPaddingBitsOfType
                          );
        uint8_t current_rhs_bit 
            = GetStartBit(
                          GetShift(rhs) - i - 1 - kPaddingBitsOfType
                         );
        uint8_t current_rhs_byte 
            = GetStartByte(
                           GetShift(rhs) - i - 1 - kPaddingBitsOfType
                          );
        uint8_t lhs_bit_value 
            = GetBitValue(
                          lhs
                          , current_lhs_byte
                          , current_lhs_bit
                         );
        uint8_t rhs_bit_value
            = GetBitValue(
                          rhs
                          , current_rhs_byte
                          , current_rhs_bit
                         );

        if (lhs_bit_value != rhs_bit_value)
            return lhs_bit_value < rhs_bit_value;
    }    

    return false;
}

bool operator>(const uint239_t& lhs, const uint239_t& rhs) {
    return rhs < lhs;
}

bool operator>=(const uint239_t& lhs, const uint239_t& rhs) {
    return !(lhs < rhs);
}

bool operator<=(const uint239_t& lhs, const uint239_t& rhs) {
    return !(lhs > rhs);
}

uint239_t FromInt(uint32_t value, uint64_t shift) {
    uint239_t number = kUint239Zero;

    uint8_t current_number_bit = GetStartBit(shift);
    uint8_t current_number_byte = GetStartByte(shift);

    const uint8_t UINT32_SIZE = 32;
    for (size_t i = 0; i < UINT32_SIZE; ++i) {
        uint8_t current_value_bit = value & 1;
        number.data[current_number_byte]
            += (1 << current_number_bit) * current_value_bit;
        value >>= 1;
        
        MoveToNextBit(current_number_bit, current_number_byte);
    }

    shift %= kShiftMaxValue + 1ll;
    SetShift(number, shift);
    
    return number;
}

uint239_t FromString(const char* str, uint64_t shift) {
    uint239_t number = kUint239Zero;

    const size_t STR_LEN = strlen(str);
    uint8_t binary_num[kSignificantBitsOfType]{};
    uint8_t buffer[kSignificantBitsOfType]{};
    uint8_t num[kSignificantBitsOfType]{};

    for (size_t i = 0; i < STR_LEN; ++i)
        num[i] = str[i] - '0';

    for (size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t has_ten = 0;
        binary_num[i] = num[STR_LEN - 1] & 1;
        for (size_t i = 0; i < STR_LEN; ++i) {
            uint8_t intermediate_value = has_ten * 10 + num[i];
            buffer[i] = (intermediate_value >> 1);
            has_ten = intermediate_value & 1;
        }

        for (size_t i = 0; i < STR_LEN; ++i) {
            num[i] = buffer[i];
        }
    }

    shift %= kShiftMaxValue + 1ll;
    SetShift(number, shift);

    uint8_t current_number_bit = GetStartBit(shift);
    uint8_t current_number_byte = GetStartByte(shift);
    for(size_t i = 0; i < kSignificantBitsOfType; ++i) {
        number.data[current_number_byte]
            += (1 << current_number_bit) * binary_num[i];

        MoveToNextBit(current_number_bit, current_number_byte);
    }

    return number;
}

uint239_t operator+(const uint239_t& lhs, const uint239_t& rhs) {
    uint64_t shift 
        = (
           GetShift(lhs) + GetShift(rhs)
          ) 
           % (kShiftMaxValue + 1ll);
    uint239_t number = FromInt(0, shift);

    uint8_t current_lhs_bit = GetStartBit(GetShift(lhs));
    uint8_t current_lhs_byte = GetStartByte(GetShift(lhs));
    uint8_t current_rhs_bit = GetStartBit(GetShift(rhs));
    uint8_t current_rhs_byte = GetStartByte(GetShift(rhs));
    uint8_t current_number_bit = GetStartBit(shift);
    uint8_t current_number_byte = GetStartByte(shift);
    uint8_t carry_bit = 0;

    for (size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t bits_sum 
            = ((lhs.data[current_lhs_byte] >> current_lhs_bit) & 1)
              + ((rhs.data[current_rhs_byte] >> current_rhs_bit) & 1)
              + carry_bit;
        uint8_t result_bit = bits_sum & 1;

        number.data[current_number_byte] 
            += (1 << current_number_bit) * result_bit;
        carry_bit = (bits_sum >> 1);

        MoveToNextBit(current_lhs_bit, current_lhs_byte);
        MoveToNextBit(current_rhs_bit, current_rhs_byte);
        MoveToNextBit(current_number_bit, current_number_byte);
    }

    return number;
}

uint239_t operator-(const uint239_t& lhs, const uint239_t& rhs) {
    if (rhs > lhs) {
        return kUint239Zero;
    }

    uint64_t shift 
        = (
           kShiftMaxValue + 1ll + GetShift(lhs) - GetShift(rhs)
          )
           % (kShiftMaxValue + 1ll);
    uint239_t number = FromInt(0, shift);
    
    uint8_t current_lhs_bit = GetStartBit(GetShift(lhs));
    uint8_t current_lhs_byte = GetStartByte(GetShift(lhs));
    uint8_t current_rhs_bit = GetStartBit(GetShift(rhs));
    uint8_t current_rhs_byte = GetStartByte(GetShift(rhs));
    uint8_t current_number_bit = GetStartBit(shift);
    uint8_t current_number_byte = GetStartByte(shift);
    uint8_t carry_bit = 0;

    for (size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t lhs_bit_value 
            = GetBitValue(
                          lhs
                          , current_lhs_byte
                          , current_lhs_bit
                         );
        uint8_t rhs_bit_value 
            = GetBitValue(
                          rhs
                          , current_rhs_byte
                          , current_rhs_bit
                         );
        uint8_t number_bit_value 
            = (2 + lhs_bit_value - rhs_bit_value - carry_bit) & 1;
        
        carry_bit 
            = (lhs_bit_value < rhs_bit_value + carry_bit);
        number.data[current_number_byte] 
            += number_bit_value << current_number_bit;

        MoveToNextBit(current_lhs_bit, current_lhs_byte);  
        MoveToNextBit(current_rhs_bit, current_rhs_byte);
        MoveToNextBit(current_number_bit, current_number_byte);
    }

    return number;
}

uint239_t operator*(const uint239_t& lhs, const uint239_t& rhs) {
    uint64_t shift = (GetShift(lhs) + GetShift(rhs)) % (kShiftMaxValue + 1ll);
    uint239_t number = FromInt(0, shift);
    
    uint8_t current_rhs_bit = GetStartBit(GetShift(rhs));
    uint8_t current_rhs_byte = GetStartByte(GetShift(rhs));

    for (size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t rhs_bit_value 
            = GetBitValue(
                          rhs
                          , current_rhs_byte
                          , current_rhs_bit
                         );
        MoveToNextBit(current_rhs_bit, current_rhs_byte);
        
        if (rhs_bit_value != 1)
            continue;

        uint239_t buffer = lhs;
        ChangeShift(buffer, 0);
        MoveNumberToLeft(buffer, i);
        number = number + buffer;
    }

    return number;
}

uint239_t operator/(const uint239_t& lhs, const uint239_t& rhs) {  
    uint64_t shift 
        = (kShiftMaxValue + 1ll + GetShift(lhs) - GetShift(rhs)) 
           % (kShiftMaxValue + 1ll);
    uint239_t number = FromInt(0, shift);

    uint239_t lhs_copy = lhs;
    uint239_t rhs_copy = rhs;
    ChangeShift(lhs_copy, 0);
    ChangeShift(rhs_copy, 0);

    uint8_t max_divider = 0;
    while (lhs_copy >= rhs_copy) {
        MoveNumberToLeft(rhs_copy, 1);
        ++max_divider;
    }
    if (max_divider > 0){
        --max_divider;
    }

    for (size_t j = 0; j <= max_divider; ++j) {
        number = number * FromInt(2, 0);
        rhs_copy = rhs;
        MoveNumberToLeft(rhs_copy, max_divider - j);

        if (rhs_copy > lhs_copy) 
            continue;

        lhs_copy = lhs_copy - rhs_copy;
        number = number + FromInt(1, 0);
    }

    return number;
}

bool operator==(const uint239_t& lhs, const uint239_t& rhs) {
    uint8_t current_lhs_bit = GetStartBit(GetShift(lhs));
    uint8_t current_lhs_byte = GetStartByte(GetShift(lhs));
    uint8_t current_rhs_bit = GetStartBit(GetShift(rhs));
    uint8_t current_rhs_byte = GetStartByte(GetShift(rhs));

    for (size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t lhs_bit_value 
            = GetBitValue(
                          lhs
                          , current_lhs_byte
                          , current_lhs_bit
                         );
        uint8_t rhs_bit_value 
            = GetBitValue(
                          rhs
                          , current_rhs_byte
                          , current_rhs_bit
                         ); 
        
        if (lhs_bit_value != rhs_bit_value) {
            return false;
        }

        MoveToNextBit(current_lhs_bit, current_lhs_byte);
        MoveToNextBit(current_rhs_bit, current_rhs_byte);
    }    

    return true;
}

bool operator!=(const uint239_t& lhs, const uint239_t& rhs) {
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& stream, const uint239_t& value) {
    uint8_t buffer[kSignificantBitsOfType]{};
    uint8_t output[kSignificantBitsOfType]{};
    buffer[0] = 1;

    uint8_t current_bit = GetStartBit(GetShift(value));
    uint8_t current_byte = GetStartByte(GetShift(value));

    for (size_t i = 0; i < kSignificantBitsOfType; ++i) { 
        if (GetBitValue(value, current_byte, current_bit) == 1) {
            uint8_t carry_digit = 0;
            for (size_t j = 0; j < kSignificantBitsOfType; ++j) {
                uint8_t digit_sum = output[j] + buffer[j] + carry_digit;
                output[j] = digit_sum % 10;
                carry_digit = digit_sum / 10;
            }
        }

        uint8_t carry_digit = 0;
        for (size_t j = 0; j < kSignificantBitsOfType; ++j) {
            uint8_t digit_sum = (buffer[j] >> 1) + carry_digit;
            buffer[j] = digit_sum % 10;
            carry_digit = digit_sum / 10;
        }

        MoveToNextBit(current_bit, current_byte);
    }

    bool is_insignificant_zero = true;
    for (size_t i = 0; i < kSignificantBitsOfType; ++i) {
        uint8_t current_digit = output[kSignificantBitsOfType - 1 - i];
        if (is_insignificant_zero) {
            if (current_digit == 0 && i != kSignificantBitsOfType - 1) {
                continue;
            }
            is_insignificant_zero = false;
        }
        stream << static_cast<int>(current_digit);
    }
    stream << ", shift = " << GetShift(value);

    return stream;
}

uint64_t GetShift(const uint239_t& value) {
    uint64_t shift = 0;
    for (size_t current_byte = 0; current_byte < kUint239Size; ++current_byte) {
        shift <<= 1;
        shift += (value.data[current_byte] >> 7);
    }

    return shift;
}
