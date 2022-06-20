#pragma once
#include <vector>
#include <cstdint>
#include <iostream>
#include <string>
#include <algorithm>
#include <intrin.h>
#include <cmath>

#pragma intrinsic(_umul128)


template <size_t n>
class BigFloat
{
	bool sign = false;
	int32_t exp = 0;
    const static inline std::vector<std::vector<BigFloat<n>>> *powers_of_10 = nullptr;
    const static inline BigFloat<n> *ten = nullptr;
    const static inline BigFloat<n> *zero = nullptr;
    const static inline BigFloat<n> *one_c = nullptr;
    const static inline uint32_t log_2 = 2585827973;
public:

    std::vector<std::uint64_t > mantisse;

    BigFloat() :
        sign(false),
        exp(0),
        mantisse()
    {
        mantisse.resize(n);
    }


    BigFloat(const double x) :
        sign(false),
        exp(0),
        mantisse()
    {
        mantisse.resize(n);

        uint64_t mant = 0;
        if (x != 0)
        {
            uint64_t u64 = *reinterpret_cast<const uint64_t*>(&x);
            mant = (u64 << 11) | (1ULL << 63);
            sign = u64 >> 63;
            exp = (static_cast<int64_t>(u64 >> 52) & 0x7FF) - 1023;
        }
        mantisse.resize(n);
        mantisse[n-1] = mant;      
    }

    
    BigFloat(std::vector<std::uint64_t > p_mantisse, bool p_sign, int32_t p_exp)
    {
        sign = p_sign;
        exp = p_exp;
        std::reverse(p_mantisse.begin(), p_mantisse.end());
        mantisse = p_mantisse;
    }


    BigFloat(const std::string& str) :
        sign(false),
        exp(0),
        mantisse()
    {
        mantisse.resize(n);

        uint64_t mant = 0;
        bool negativesign = false;
        size_t i = 0;

        if(str.empty())
            throw std::invalid_argument("string could not be converted to BigFloat");


        if (str[i] == '-')
            negativesign = true, i++;
        else if (str[i] == '+')
            i++;

        const BigFloat<n> BF10 = BigFloat(10.0);
        size_t numpredigits = 0;
        while (i < str.size() && isdigit(str[i])) 
        {
            *this = *this * BF10;
            *this = *this + BigFloat(static_cast<double>(str[i] - '0'));
            i++;
            numpredigits++;
        }

        int32_t numpostdigits = 0;
        int32_t exp_after_e = 0;
        if(i < str.size() && str[i] == '.')
        {
            i++;

            while (i < str.size() && isdigit(str[i]))
            {
                numpostdigits++;
                *this = *this * BF10;
                *this = *this + BigFloat(static_cast<double>(str[i] - '0'));
                i++;
            }
        }

        if (i < str.size() && (str[i] == 'e' || str[i] == 'E'))
        {
            i++;
            size_t numExpDigits = 0;
            if (i < str.size()) 
            {
                size_t indexBehindNumber = 0;
                try
                {
                    exp_after_e = stol(str.substr(i, str.size()), &indexBehindNumber);
                }
                catch (std::out_of_range) 
                {
                    throw std::invalid_argument("exponent too large for BigFloat");
                }
                catch (...) 
                {
                    throw;
                }
                numExpDigits = str.size() - i;
                i += indexBehindNumber;
            }
            else
                throw std::invalid_argument("string could not be converted to BigFloat");
        }

        if (powers_of_10 == nullptr)
            create_table();

        const BigFloat<n> correctionFactor = pow_10(abs(numpostdigits - exp_after_e));
        if ((numpostdigits - exp_after_e) < 0)
            *this = *this * correctionFactor;
        else
            *this = *this / correctionFactor;
        if (i != str.size())
            throw std::invalid_argument("string could not be converted to BigFloat");
        
    }

    int32_t getExp() { return exp; }

    static void sqrt(BigFloat& res)
    {
        if (res.sign)
            return;
        BigFloat<n> num = res;
        int itterations = (int)std::sqrt(n * 22.0);
        BigFloat<n> oneOverTwo = 0.5;
        BigFloat<n> numOverRes = 0.0;
        BigFloat<n> temp;
        for (int i = 0; i <= itterations; i++)
        {
            real_divide(numOverRes, num, res);
            temp = res;
            real_add(res, temp, numOverRes);
            temp = res;
            multiply(res, oneOverTwo, temp);
            std::cout << "\n" << (100*i / itterations);
        }
    }

    void print()
    {
        int32_t exp_base10 = static_cast<int64_t>(abs(exp)) * static_cast<int64_t>(log_2) >> 33;


        if (powers_of_10 == nullptr)
            create_table();

        if (ten == nullptr)
            create_constants();

        BigFloat correction_factor = pow_10(exp_base10);
        BigFloat temp = *this;
        BigFloat rounder(1);
        int counter = static_cast<int64_t>(log_2) * (64 * n) >> 33;
        rounder.exp = -static_cast<int32_t>((n*64-1));

        if (exp < 0)
            multiply(temp, *this, correction_factor);
        else
        {
            real_divide(temp, *this, correction_factor);
        }

        bool extra_step = false;

        if (compare_mant(temp, *ten) == 1)
        {
            temp = temp / *ten; 
            extra_step = true;
        }
        else if (compare_mant(temp, *one_c) == 2)
        {
            multiply(temp, temp, *ten);
            extra_step = true;
        }
            
        if (!sign)
            temp = temp + rounder;
        else
            temp = temp - rounder;

        unsigned char digit = 0;
        std::string decimal = "";

        digit = temp.mantisse[n - 1] >> (63 - temp.exp);
        //digit = temp.mantisse[n - 1] >> (63 - temp.exp-1);

        decimal += '0' + digit;

        shift_inner_l(temp, temp.exp+1);

        BigFloat left_shift_3 = temp;
        BigFloat left_shift_1 = temp;
        

        while (counter > 0)
        {
            digit = temp.multiply_10(left_shift_1, left_shift_3);
            decimal += '0' + digit;
            left_shift_3 = temp;
            left_shift_1 = temp;
            counter--;
        }
        if (!extra_step)
            decimal.insert(exp_base10 + 1, ".");
        else
            decimal.insert(exp_base10 + 2, ".");
        //std::reverse(decimal.begin(), decimal.end());
        std::cout << decimal;

    }


private:

    template <size_t>
    friend class BigFloat;

    unsigned char multiply_10(BigFloat& left_shift_1,BigFloat& left_shift_3)
    {
        shift_inner_l(left_shift_1, 1);
        shift_inner_l(left_shift_3, 3);
        unsigned char digit = (mantisse[n - 1] >> 63) + (mantisse[n - 1] >> 61);

        unsigned char carry = 0;
        for (size_t i = 0; i < n; i++)
        {
            carry = _addcarry_u64 (carry, left_shift_1.mantisse[i], left_shift_3.mantisse[i], &mantisse[i]);
        }

        digit += carry;
        return digit;

    }

    static BigFloat pow_10(int32_t e)
    {
        BigFloat res(1);
        int32_t temp;
        size_t i = 0;
        while (e != 0)
        {
            temp = e & 0xF;
            res = res * (*powers_of_10)[i][temp];
            e >>= 4;
            i++;             
        }
        return res;
    }

    static void create_table()
    {
        BigFloat<n + 1> temp = 10;
        BigFloat<n + 1> base = 10;
        BigFloat<n> small;
        BigFloat<n> one(1);
        auto p  = new std::vector<std::vector<BigFloat<n>>> (8, std::vector<BigFloat<n>>(16));

        for (size_t i = 0; i < 8; i++)
        {
            (*p)[i][0] = one;
            for (size_t j = 1; j < 16; j++)
            {
                small.copy_cutoff(temp);
                (*p)[i][j] = small;
                temp = temp * base;
            }
            base = temp;
        }
        powers_of_10 = p;
    }

    static void create_constants()
    {
        ten = new BigFloat(10.0);
        zero = new BigFloat(0);
        one_c = new BigFloat(1);
    }

    static void shift(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs, const bool lhs_greater)
    {
        const size_t shift_val = std::abs(lhs.exp - rhs.exp);

        if (lhs_greater)
        {
            res = rhs;
        }
        else
        {
            res = lhs;
        }
        shift_inner_r(res, shift_val);
    }

    void copy_cutoff(const BigFloat<n+1>& val)
    {
        for (size_t i = n-1; i > 0; i--)
        {
            mantisse[i] = val.mantisse[i+1];
        }
        exp = val.exp;
        return;
    }

    static void shift_inner_r(BigFloat& res, const size_t shift_val)
    {
        const size_t source = shift_val / 64;
        const size_t first_shift = (shift_val % 64);
        const size_t follow_shift = 64 - first_shift;

        if (source < n)
        {
            if (first_shift == 0)
            {
                for (size_t i = 0; i < (n-source); i++)
                {
                    res.mantisse[i] = res.mantisse[i + source];
                }
                
            }
            else
            {
                for (size_t i = 0; i < n - 1 - (source); i++)
                {
                    res.mantisse[i] = res.mantisse[i + source] >> first_shift;
                    res.mantisse[i] |= res.mantisse[i + source + 1] << follow_shift;
                }
                res.mantisse[n - 1 - source] = res.mantisse[n - 1 - source] >> first_shift;
                if (source > 0)
                    res.mantisse[n - 1 - source] = res.mantisse[n - 1] >> first_shift;
            }
            std::fill(res.mantisse.begin() + n - source, res.mantisse.end(), 0);
            
        }
        else
            std::fill(res.mantisse.begin(), res.mantisse.end(), 0);
    }
    static void shift_inner_l(BigFloat& res, const size_t shift_val)
    {
        const size_t source = shift_val / 64;
        const size_t first_shift = (shift_val % 64);
        const size_t follow_shift = 64 - first_shift;

        if (source < n)
        {
            if (first_shift == 0)
            {
                for (size_t i = n-1; i > source; i--)
                {
                    res.mantisse[i] = res.mantisse[i - source] << first_shift;
                }

            }
            else
            {
                for (size_t i = n - 1; i > source; i--)
                {
                    res.mantisse[i] = res.mantisse[i - source] << first_shift;
                    res.mantisse[i] |= res.mantisse[i - source - 1] >> follow_shift;
                }
            }
            res.mantisse[source] = res.mantisse[source] << first_shift;
            if (source > 0)
                std::fill(res.mantisse.begin(), res.mantisse.begin()+source-1, 0);

        }
        else
            std::fill(res.mantisse.begin(), res.mantisse.end(), 0);
    }

    static void shift_format(BigFloat& res)
    {
        int shift_val = 0;
        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = 63; j >= 0; j--)
            {
                if (res.mantisse[i] < (1ULL << j))
                    shift_val++;
                else
                {
                    i = -1;
                    break;
                }
            }
        }
        shift_inner_l(res, shift_val);
        res.exp -= shift_val;
    }

    static void negate(BigFloat& res)
    {
        bool carryBit = false;
        for (size_t i = 0; i < n; i++)
        {
            res.mantisse[i] = 0 - res.mantisse[i] - carryBit;
            if (res.mantisse[i] != 0)
                carryBit = 1;
            else
                carryBit = 0;
        }
    }

    static size_t compare_mant(const BigFloat& lhs, const BigFloat& rhs)
    {
        // 0 = equal
        // 1 = lhs is greater
        // 2 = rhs is greater

        if (lhs.exp > rhs.exp)
            return 1;
        else if (lhs.exp < rhs.exp)
            return 2;
        else
        {
            for (int i = n - 1; i >= 0; i--)
            {
                if (lhs.mantisse[i] > rhs.mantisse[i])
                    return 1;
                else if (lhs.mantisse[i] < rhs.mantisse[i])
                    return 2;
            }
            return 0;
        }
    }
    static size_t mant_compare(const BigFloat& lhs, const BigFloat& rhs)
    {
        for (int i = n - 1; i >= 0; i--)
        {
            if (lhs.mantisse[i] > rhs.mantisse[i])
                return 1;
            else if (lhs.mantisse[i] < rhs.mantisse[i])
                return 2;
        }
        return 0;
    }

    static void mant_add(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs)
    {
        unsigned char carryBit = 0;
        res.sign = lhs.sign;
        if (lhs.exp > rhs.exp)
        {
            shift(res, lhs, rhs, true);
            res.exp = lhs.exp;

            for (size_t i = 0; i < n; i++)
            {
                carryBit = _addcarry_u64(carryBit, res.mantisse[i], lhs.mantisse[i], &res.mantisse[i]);
            }
        }
        else
        {
            shift(res, lhs, rhs, false);
            res.exp = rhs.exp;

            for (size_t i = 0; i < n; i++)
            {
                carryBit = _addcarry_u64(carryBit, res.mantisse[i], rhs.mantisse[i], &res.mantisse[i]);
            }
        }

        if (carryBit)
        {
            shift_inner_r(res, 1);
            res.mantisse[n - 1] |= 1ULL << 63;
            res.exp++;
        }
    }

    void mant_add_u64(uint64_t val, size_t index)
    {
        bool carryBit = false;
        do 
        {
            uint64_t temp = mantisse[index];
            mantisse[index] += val + carryBit;
            val = 0;
            if ((mantisse[index] < temp && !carryBit) || (mantisse[index] == temp && carryBit))
            {
                carryBit = 1;
                index++;
            }
            else
                carryBit = 0;
        } 
        while (carryBit);

    }

    static void multiply(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs)
    {
        if (rhs.mantisse[n - 1] >> 63 && lhs.mantisse[n - 1] >> 63)
        {
            uint64_t hi = 0; uint64_t lo = 0;
            BigFloat<2 * n> res_temp = 0;
            bool carryBit = false;

            for (size_t l = 0; l < n; l++)
            {
                if (lhs.mantisse[l] != 0)
                {
                    for (size_t r = 0; r < n; r++)
                    {
                        lo = _umul128(lhs.mantisse[l], rhs.mantisse[r], &hi);
                        res_temp.mant_add_u64(lo, l + r);
                        res_temp.mant_add_u64(hi, l + r + 1);
                    }
                }
            }
            res.exp = lhs.exp + rhs.exp;
            for (int i = n - 1; i >= 0; i--)
            {
                res.mantisse[i] = res_temp.mantisse[n + i];
            }
            if (res_temp.mantisse[2 * n - 1] >> 63)
                res.exp += 1;
            else
            {
                shift_inner_l(res, 1);
                res.mantisse[0] += res_temp.mantisse[n - 1] >> 63;
            }
            res.sign = lhs.sign ^ rhs.sign;
        }   
        else if (res.mantisse[n - 1] >> 63)
        {
            res.sign = false;
            std::fill(res.mantisse.begin(), res.mantisse.end(), 0);
        }
        
    }

    bool is_zero() const
    {
        for (size_t i = 0; i < n; i++)
        {
            if (mantisse[i] != 0)
                return false;
        }
        return true;
    }

    void set_bit(size_t index,const bool bit)
    {
        size_t index_bit = index % 64;
        index = index / 64;
        if (bit)
            mantisse[index] |= 1ULL << index_bit;
        else
            mantisse[index] &= ~(1ULL << index_bit);
    }

    static bool div_subtract(BigFloat& temp, const BigFloat& rhs, bool& overflow_bit, bool& end)
    {
        if (!overflow_bit)
        {
            const size_t cmp = mant_compare(temp, rhs);

            if (cmp == 1)    // temp > rhs
            {
                unsigned char carry = 0;
                for (size_t i = 0; i < n; i++)
                {
                    carry = _subborrow_u64(carry, temp.mantisse[i], rhs.mantisse[i], &temp.mantisse[i]);
                }
                shift_inner_l(temp, 1);
                return false;
            }
            else if (cmp == 2)  // rhs > temp
            {
                if (temp.mantisse[n - 1] >> 63)       // if first bit = 1
                    overflow_bit = true;
                shift_inner_l(temp, 1);
                return true;
            }
            else
            {
                end = true;
                return false;
            }    
        }
        else
        {
            unsigned char carry = 0;
            for (size_t i = 0; i < n; i++)
            {
                carry = _subborrow_u64(carry, temp.mantisse[i], rhs.mantisse[i], &temp.mantisse[i]);
            }
            if (!(temp.mantisse[n - 1] >> 63))       // if first bit = 0
                overflow_bit = false;
            shift_inner_l(temp, 1);
            return false;
        }

        //for (size_t i = 0; i < n; i++)
        //{
        //    res.mantisse[i] = temp.mantisse[i] - rhs.mantisse[i] - carryBit;
        //    carryBit = 0;
        //    if (res.mantisse[i] > temp.mantisse[i] || ((res.mantisse[i] == temp.mantisse[i]) && rhs.mantisse[i] != 0))               //|| ((res.mantisse[i] == lhs.mantisse[i]) && temp != 0)
        //    {
        //        carryBit = 1;     
        //    }
        //    end |= res.mantisse[i];
        //}

        //if (!overflow_bit)
        //{
        //   
        //    if (carryBit)
        //    {
        //        if (temp.mantisse[n - 1] >> 63)       // if first bit = 1
        //            overflow_bit = true;

        //        shift_inner_l(temp, 1);
        //        return true;
        //    }
        //    else
        //    {
        //        temp = res;
        //        shift_inner_l(temp, 1);
        //        return false;
        //    }

        //    
        //}
        //else
        //{
        //    temp = res;
        //    if (!(temp.mantisse[n - 1] >> 63))       // if first bit != 1
        //        overflow_bit = false;      
        //    shift_inner_l(temp, 1);
        //    return false;

//}


    }

    static void divide(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs)
    {
        BigFloat<n> temp = lhs;
        int index = n * 64 - 1;
        bool preceding_overflow = false;
        bool end = false;

        do
        {
            if (div_subtract(temp, rhs, preceding_overflow, end))
                res.set_bit(index, 0);
            else
                res.set_bit(index, 1);
            index--;
        } while (!end && index >= 0);
    }

    static void mant_subtract(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs)
    {
        bool carryBit = 0;
        const size_t cmp = compare_mant(lhs, rhs);

        if (cmp == 1) // left > right
        {
            shift(res, lhs, rhs, true);
            res.exp = lhs.exp;
            for (size_t i = 0; i < n; i++)
            {
                uint64_t temp = res.mantisse[i];
                res.mantisse[i] = lhs.mantisse[i] - temp - carryBit;
                carryBit = 0;
                if ((res.mantisse[i] > lhs.mantisse[i]) || ((res.mantisse[i] == lhs.mantisse[i]) && temp != 0))
                {
                    carryBit = 1;
                }

            }
        }
        else if (cmp == 2) // right > left
        {
            shift(res, lhs, rhs, false);
            res.exp = rhs.exp;
            for (size_t i = 0; i < n; i++)
            {
                uint64_t temp = res.mantisse[i];
                res.mantisse[i] = temp - rhs.mantisse[i] - carryBit;
                carryBit = 0;
                if ((res.mantisse[i] > temp) || ((res.mantisse[i] == temp) && rhs.mantisse[i] != 0))
                {
                    carryBit = 1;
                }
            }
        }
        else
            return;

        if (carryBit)
        {
            negate(res);
            res.sign = !(lhs.sign);
        }
        else
            res.sign = lhs.sign;

        shift_format(res);
    }

    static void real_divide(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs)
    {
        if (!rhs.is_zero())
        {
            if (!lhs.is_zero())
            {
                divide(res, lhs, rhs);
                res.sign = lhs.sign ^ rhs.sign;
                res.exp = lhs.exp - rhs.exp;
                shift_format(res);
            }
            else
                res = lhs;
        }
    }

    static void real_add(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs)
    {
        if ((lhs.mantisse[n - 1] >> 63) == 0)      //if left = 0
            res = rhs;
        else if ((rhs.mantisse[n - 1] >> 63) == 0)   //if right = 0
            res = lhs;


        if (lhs.sign == rhs.sign)
        {
            mant_add(res, lhs, rhs);
        }
        else
            mant_subtract(res, lhs, rhs);
    }

    static void real_subtract(BigFloat& res, const BigFloat& lhs, const BigFloat& rhs)
    {
        if (lhs.sign != rhs.sign)
        {
            mant_add(res, lhs, rhs);
            res.sign = lhs.sign;
        }
        else
        {
            mant_subtract(res, lhs, rhs);
        }
    }


    friend BigFloat<n> operator/(const BigFloat& lhs, const BigFloat& rhs)
    {
        BigFloat<n> res(0);
        if (!rhs.is_zero())
        {
            if (!lhs.is_zero())
            {
                divide(res, lhs, rhs);
                res.sign = lhs.sign ^ rhs.sign;
                res.exp = lhs.exp - rhs.exp;
                shift_format(res);
                return res;
            }
            else
                return res;
        }
        //else
            //exception
        return res;
    }

    friend BigFloat<n> operator+(const BigFloat& lhs, const BigFloat& rhs)
    {
        BigFloat<n> res(0);
        if ((lhs.mantisse[n - 1] >> 63) == 0)      //if left = 0
            return rhs;
        else if ((rhs.mantisse[n - 1] >> 63) == 0)   //if right = 0
            return lhs;


        if (lhs.sign == rhs.sign)
        {
            mant_add(res, lhs, rhs);
        }
        else
            mant_subtract(res, lhs, rhs);

        return res;
    }

    friend BigFloat<n> operator-(const BigFloat& lhs, const BigFloat& rhs)
    {
        BigFloat<n> res(0);

        if (lhs.sign != rhs.sign)
        {
            mant_add(res, lhs, rhs);
            res.sign = lhs.sign;
        }
        else
        {
            mant_subtract(res, lhs, rhs);
        }
        return res;

    }



    friend BigFloat<n> operator*(const BigFloat& lhs, const BigFloat& rhs)
    {
        BigFloat<n> res(0);
        multiply(res, lhs, rhs);
        return res;
    }

};


