#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#define SIZE 512
int main() {
	unsigned char test[SIZE] = "";

	time_t timep;
	int i = 0;

	struct tm *ptr;
	timep = time(NULL);
	ptr = localtime(&timep);
	strftime(test,50, "%G%m%d-%H%M%S\n", ptr);
	printf("%s", test);
//now time: 20171103172430 输出字符格式(atoi转换为数字)
	return 0;
}
