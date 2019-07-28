#ifndef __STRUTIL_HH
#define __STRUTIL_HH

#include <string>
#include <vector>
#include <sstream>

namespace StrUtil
{
	template <class T> std::string sT(const T &t)
	{
		std::ostringstream ss;
		ss << t;
		return ss.str();
	};

	void replace_substring(std::string& _s, const std::string& old_s, const std::string& new_s);
	std::string replace_substring_copy(std::string& _s, const std::string& old_s, const std::string& new_s);

	void replace_substring(std::string& _s, const char* old_s, const char* new_s);
	std::string replace_substring_copy(std::string& _s, const char* old_s, const char* new_s);

	double str_to_double(const char *str, double defret = 0.0);
	long str_to_long(const char *str, long defret = 0);
	int str_to_int(const char *str, int defret = 0);

	std::vector<std::string> split_to_vector(const char *cs, const char *csdelim);

	class Args
	{
	public:
		Args();
		Args(const char *s);
		Args(const std::string &s);
		const char * get(const char *swhat, const char *sdefault = "");
		double getD(const char *swhat, double default_value = 0);
		int getI(const char *swhat, int default_value = 0);
		~Args();
	private:
		std::string _s;
	};
};

#endif //  __STRUTIL_HH
