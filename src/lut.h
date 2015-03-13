#ifndef __LUT_H
#define __LUT_H

const unsigned char base_to_number[] = {
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
99,
0,
99,
1,
99,
99,
99,
2,
99,
99,
99,
99,
99,
99,
4,
99,
99,
99,
99,
99,
3,
99,
99,
99,
99,
99,
99
};

// for (unsigned char i = 0; i <= 90; i ++)
// {
// 	if (i == 'A') cout << "0," << endl;
// 	else if (i == 'C') cout << "1," << endl;
// 	else if (i == 'G') cout << "'G'," << endl;
// 	else if (i == 'T') cout << "'T'," << endl;
// 	else if (i == 'N') cout << "'N'," << endl;
// 	else cout << "99," << endl;
// }

const char number_to_basetriplet[125][3] = {
{'A','A','A'},
{'A','A','C'},
{'A','A','G'},
{'A','A','T'},
{'A','A','N'},
{'A','C','A'},
{'A','C','C'},
{'A','C','G'},
{'A','C','T'},
{'A','C','N'},
{'A','G','A'},
{'A','G','C'},
{'A','G','G'},
{'A','G','T'},
{'A','G','N'},
{'A','T','A'},
{'A','T','C'},
{'A','T','G'},
{'A','T','T'},
{'A','T','N'},
{'A','N','A'},
{'A','N','C'},
{'A','N','G'},
{'A','N','T'},
{'A','N','N'},
{'C','A','A'},
{'C','A','C'},
{'C','A','G'},
{'C','A','T'},
{'C','A','N'},
{'C','C','A'},
{'C','C','C'},
{'C','C','G'},
{'C','C','T'},
{'C','C','N'},
{'C','G','A'},
{'C','G','C'},
{'C','G','G'},
{'C','G','T'},
{'C','G','N'},
{'C','T','A'},
{'C','T','C'},
{'C','T','G'},
{'C','T','T'},
{'C','T','N'},
{'C','N','A'},
{'C','N','C'},
{'C','N','G'},
{'C','N','T'},
{'C','N','N'},
{'G','A','A'},
{'G','A','C'},
{'G','A','G'},
{'G','A','T'},
{'G','A','N'},
{'G','C','A'},
{'G','C','C'},
{'G','C','G'},
{'G','C','T'},
{'G','C','N'},
{'G','G','A'},
{'G','G','C'},
{'G','G','G'},
{'G','G','T'},
{'G','G','N'},
{'G','T','A'},
{'G','T','C'},
{'G','T','G'},
{'G','T','T'},
{'G','T','N'},
{'G','N','A'},
{'G','N','C'},
{'G','N','G'},
{'G','N','T'},
{'G','N','N'},
{'T','A','A'},
{'T','A','C'},
{'T','A','G'},
{'T','A','T'},
{'T','A','N'},
{'T','C','A'},
{'T','C','C'},
{'T','C','G'},
{'T','C','T'},
{'T','C','N'},
{'T','G','A'},
{'T','G','C'},
{'T','G','G'},
{'T','G','T'},
{'T','G','N'},
{'T','T','A'},
{'T','T','C'},
{'T','T','G'},
{'T','T','T'},
{'T','T','N'},
{'T','N','A'},
{'T','N','C'},
{'T','N','G'},
{'T','N','T'},
{'T','N','N'},
{'N','A','A'},
{'N','A','C'},
{'N','A','G'},
{'N','A','T'},
{'N','A','N'},
{'N','C','A'},
{'N','C','C'},
{'N','C','G'},
{'N','C','T'},
{'N','C','N'},
{'N','G','A'},
{'N','G','C'},
{'N','G','G'},
{'N','G','T'},
{'N','G','N'},
{'N','T','A'},
{'N','T','C'},
{'N','T','G'},
{'N','T','T'},
{'N','T','N'},
{'N','N','A'},
{'N','N','C'},
{'N','N','G'},
{'N','N','T'},
{'N','N','N'}
};

// for (int i = 0; i < 5; i++)
// 	for (int j = 0; j < 5; j++)
// 		for (int k = 0; k < 5; k++)
// 			cout << "{" << i << "," << j << "," << k << "}," << endl;


#endif // __LUT_H