#include <Rcpp.h>
#include <bitset>
#include <vector>
#include <string>
using namespace Rcpp;


class Bit_set_class
{
public:
	size_t element_num;
	std::vector<std::bitset<8>> bit_list;
	Bit_set_class(size_t eltNum);
	size_t count_interception(const Bit_set_class& another);
	bool contain(const Bit_set_class& another);
	size_t count();
	void set_bit(std::vector<size_t>& index);
	void set_bit(size_t i);
	std::string to_string();
};

