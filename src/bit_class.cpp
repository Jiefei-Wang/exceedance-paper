#include "bit_class.h"
#include <cmath>

using std::bitset;
using std::string;
using std::vector;

Bit_set_class::Bit_set_class(size_t element_num) : element_num(element_num)
{
	int true_element = ceil(element_num / 8.0);
	bit_list.reserve(true_element);
	for (int i = 0; i < true_element; i++)
	{
		bit_list.emplace_back();
	}
}

size_t Bit_set_class::count_interception(const Bit_set_class &another)
{
	if (another.element_num != element_num)
	{
		Rf_error("The element numbers are not equal");
	}
	size_t result = 0;
	for (size_t i = 0; i < bit_list.size(); i++)
	{
		const bitset<8> &left = bit_list[i];
		const bitset<8> &right = another.bit_list[i];
		result += (left & right).count();
	}
	return result;
}

bool Bit_set_class::contain(const Bit_set_class &another)
{
	for (size_t i = 0; i < bit_list.size(); i++)
	{
		const bitset<8> &parent = bit_list[i];
		const bitset<8> &child = another.bit_list[i];
		bool result = ((~parent) & child).count() == 0;
		if (!result)
			return false;
	}
	return true;
}

void Bit_set_class::set_bit(std::vector<size_t> &index)
{
	for (auto &i : bit_list)
	{
		i.reset();
	}
	for (const auto i : index)
	{
		set_bit(i);
	}
}

void Bit_set_class::set_bit(size_t i)
{
	if (i >= element_num)
			return;
	size_t base_index = i / 8L;
	size_t off = i % 8L;
	bit_list[base_index].set(off);
}

string Bit_set_class::to_string()
{
	string str;
	for (const auto &i : bit_list)
	{
		//Rprintf("%d:%s\n", i, bit_list[i]->to_string().c_str());
		str.append(i.to_string());
	}
	reverse(str.begin(), str.end());
	return str;
}

size_t Bit_set_class::count()
{
	size_t result = 0;
	for (const auto &i : bit_list)
	{
		result += i.count();
	}
	return result;
}
