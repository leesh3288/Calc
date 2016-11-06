#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstring>
#include <stack>


std::stack<long double> OutSt;
std::stack<char> OpSt;

char equ[10000];

void input();
void process();
void output();

int main(void)
{
	input();
	process();
	output();
}

void input()
{
	printf("Input Equation\n >> ");
	scanf("%[^\n]", equ);

	return;
}

void process()
{
	int ite = 0, eqlen = strlen(equ);

	while (ite < eqlen)
	{
		if (equ[ite] == ' ')
			ite++;
		else if (equ[ite] >= '0' && equ[ite] <= '9')
		{
			long double inp = 0;
			int nite;
			for (nite = ite; equ[nite] >= '0'&&equ[nite] <= '9'; nite++)
			{
				inp *= 10;
				inp += (equ[nite] - '0');
			}

			if (equ[nite] == '.')
			{
			}
		}
		else if (equ[ite] >= 'a'&&equ[ite] <= 'z')
		{
		}
		else if (equ[ite] == ',')
		{
		}
		else if (equ[ite] == '(')
		{
		}
		else if (equ[ite] == ')')
		{
		}
		else // operator goes here
		{
		}
	}

	return;
}

void output()
{


	return;
}

// with functions with variable argument, need to separate somehow with indicators

/*
Operator/Function Hash Values:
1 = LParen
2 = Fact
3 = Exp
4 = Mult
5 = Div
6 = Add
7 = Sub
8~20 = Reserved
(1>2>3>4, 5>6, 7)

21 = max
22 = min

23 = sin
24 = cos
25 = tan
26 = csc
27 = sec
28 = cot

29 = asin
30 = acos
31 = atan
32 = acsc
33 = asec
34 = acot

35 = sinh
36 = cosh
37 = tanh
38 = csch
39 = sech
40 = coth

41 = asinh
42 = acosh
43 = atanh
44 = acsch
45 = asech
46 = acoth

47 = log
*/

/*
Mathematical constants are used with a preceding '#'
Lists: #pi, #e, #phi
*/