#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstring>
#include <cmath>
#include <stack>
#include <algorithm>

#define MAXL_CONST 10
#define MAXL_FUNC 10

#define CONST_E 2.71828182845904523536L
#define CONST_PI 3.14159265358979323846L
#define CONST_PHI 1.61803398874989484820L

#define EPS 1E-10 /// TODO: Testing out epsilon values for proper calc / Is is proper to use such values?

std::stack<long double> OutSt; // originally an output QUEUE, implemented as a stack for easy calculation (however, calculations must be done in reverse order)
std::stack<char> OpSt;

char equ[100000];

void input();
void process();
void output();
void OFProc(char);
int OPrecComp(char, char);
char FHash(char[]);

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
		if (equ[ite] == ' ') // ignore spacing
			ite++;
		else if (equ[ite] >= '0' && equ[ite] <= '9') // number
		{
			long double inp = 0;
			for (; equ[ite] >= '0'&&equ[ite] <= '9'; ite++)
			{
				inp *= 10;
				inp += (equ[ite] - '0');
			}

			if (equ[ite] == '.')
			{
				long double dec = 1;
				for (ite++; equ[ite] >= '0' && equ[ite] <= '9'; ite++)
				{
					dec /= 10;
					inp += (equ[ite] - '0')*dec;
				}
			}

			OutSt.push(inp);
		}
		else if (equ[ite] == '#') // mathematical constant (=number input)
		{
			char cons[MAXL_CONST + 1];
			int cct;

			for (ite++, cct = 0; (equ[ite] >= '0'&&equ[ite] <= '9') || (equ[ite] >= 'a'&&equ[ite] <= 'z'); ite++, cct++)
				cons[cct] = equ[ite];
			cons[cct] = '\0';

			if (strcmp(cons, "pi") == 0)
				OutSt.push(CONST_PI);
			else if (strcmp(cons, "e") == 0)
				OutSt.push(CONST_E);
			else if (strcmp(cons, "phi") == 0)
				OutSt.push(CONST_PHI);
		}
		else if ((equ[ite] >= 'a'&&equ[ite] <= 'z') || (equ[ite] >= 'A'&&equ[ite] <= 'Z') || equ[ite] == '_') // function
		{
			char funcs[MAXL_FUNC + 1];
			int cct;

			for (cct = 0; (equ[ite] >= '0'&&equ[ite] <= '9') || (equ[ite] >= 'a'&&equ[ite] <= 'z') || (equ[ite] >= 'A'&&equ[ite] <= 'Z') || equ[ite] == '_'; ite++, cct++)
				funcs[cct] = equ[ite];
			funcs[cct] = '\0';

			OpSt.push(FHash(funcs));
		}
		else if (equ[ite] == ',') // comma (parameter separator)
		{
			while (!OpSt.empty() && OpSt.top() != 1)
			{
				OFProc(OpSt.top());
				OpSt.pop();
			}

			if (OpSt.empty())
			{
				printf("Error: Left parenthesis not encountered, parsing failed.\n");
				return;
			}

			ite++;
		}
		else if (equ[ite] == '(') // left paren
		{
			OpSt.push(1);

			ite++;
		}
		else if (equ[ite] == ')') // right paren
		{
			while (!OpSt.empty() && OpSt.top() != 1)
			{
				OFProc(OpSt.top());
				OpSt.pop();
			}

			if (OpSt.empty())
			{
				printf("Error: Left parenthesis not encountered, parsing failed.\n");
				return;
			}

			OpSt.pop(); // pop out LParen

			if (!OpSt.empty() && OpSt.top() >= 21)
			{
				OFProc(OpSt.top()); // process function token
				OpSt.pop();
			}

			ite++;
		}
		else // operator (excludes LParen)
		{
			char opn;

			if (equ[ite] == '!') // process factorials immediately w/o pushing at OpSt
			{
				OFProc(2);

				ite++;
				continue;
			}
			else if (equ[ite] == '^') // Right-associative operators (only ^)
			{
				opn = 3;

				while (!OpSt.empty() && OpSt.top() >= 3 && OpSt.top() <= 7 && OPrecComp(opn, OpSt.top()) < 0)
				{
					OFProc(OpSt.top());
					OpSt.pop();
				}
			}
			else // Left-associative operators
			{
				if (equ[ite] == '*')
					opn = 4;
				else if (equ[ite] == '/')
					opn = 5;
				else if (equ[ite] == '+')
					opn = 6;
				else if (equ[ite] == '-')
					opn = 7;

				while (!OpSt.empty() && OpSt.top() >= 3 && OpSt.top() <= 7 && OPrecComp(opn, OpSt.top()) <= 0)
				{
					OFProc(OpSt.top());
					OpSt.pop();
				}
			}

			OpSt.push(opn);

			ite++;
		}
	}

	while (!OpSt.empty())
	{
		OFProc(OpSt.top());
		OpSt.pop();
	}

	return;
}

void output()
{
	if (OutSt.size() != 1)
	{
		printf("Error: Size of output stack != 1\n");
		return;
	}

	printf("Result: %.15lf\n", OutSt.top());

	return;
}

void OFProc(char ofnum) // operator/function processor
{
	// operator
	if (ofnum == 2) /// TODO: Floating point precision error, change to integer calc
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		int in = (int)(n + 0.5);

		n = 1;

		for (int i = 2; i <= in; i++)
			n *= i;

		OutSt.push(round(n));
	}
	else if (ofnum == 3)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(pow(n1, n2));
	}
	else if (ofnum == 4)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1*n2);
	}
	else if (ofnum == 5)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1 / n2);
	}
	else if (ofnum == 6)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1 + n2);
	}
	else if (ofnum == 7)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1 - n2);
	}

	// function
	else if (ofnum == 21)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(std::max(n1, n2));
	}
	else if (ofnum == 22)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(std::min(n1, n2));
	}
	else if (ofnum == 23)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(sin(n));
	}
	else if (ofnum == 24)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(cos(n));
	}
	else if (ofnum == 25)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(tan(n));
	}
	else if (ofnum == 26)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / sin(n));
	}
	else if (ofnum == 27)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / cos(n));
	}
	else if (ofnum == 28)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / tan(n));
	}
	else if (ofnum == 29)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asin(n));
	}
	else if (ofnum == 30)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acos(n));
	}
	else if (ofnum == 31)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atan(n));
	}
	else if (ofnum == 32)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asin(1 / n));
	}
	else if (ofnum == 33)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acos(1 / n));
	}
	else if (ofnum == 34)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atan(1 / n));
	}
	else if (ofnum == 35)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(sinh(n));
	}
	else if (ofnum == 36)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(cosh(n));
	}
	else if (ofnum == 37)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(tanh(n));
	}
	else if (ofnum == 38)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / sinh(n));
	}
	else if (ofnum == 39)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / cosh(n));
	}
	else if (ofnum == 40)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / tanh(n));
	}
	else if (ofnum == 41)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asinh(n));
	}
	else if (ofnum == 42)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acosh(n));
	}
	else if (ofnum == 43)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atanh(n));
	}
	else if (ofnum == 44)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asinh(1 / n));
	}
	else if (ofnum == 45)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acosh(1 / n));
	}
	else if (ofnum == 46)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atanh(1 / n));
	}
	else if (ofnum == 47)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(log(n2) / log(n1));
	}
	else if (ofnum == 48)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(erf(n));
	}
	else if (ofnum == 49)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(erfc(n));
	}
	else if (ofnum == 50) /// TODO: Floating point precision error, change to integer calc
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(round(tgamma(n1 + 1) / tgamma(n1 - n2 + 1)));
	}
	else if (ofnum == 51) /// TODO: Floating point precision error, change to integer calc
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(round(tgamma(n1 + 1) / (tgamma(n1 - n2 + 1) * tgamma(n2 + 1))));
	}
	else if (ofnum == 52)
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		int in1 = (int)(n1 + 0.5), in2 = (int)(n2 + 0.5);
		long long int res = 0;

		in1++; in2++;
		int *grid = new int[in1*in1]; // access by grid[in1*y+x]

		for (int i = 0; i < in1; i++)
			for (int j = 0; j < in1; j++)
				grid[in1*j + i] = 0;

		if (in1 < in2)
			in2 = 1; // forcing result to zero
		grid[in1*(in2 - 1) + in1 - 1] = 1;

		for (int i = in1 - 1; i >= 1; i--) // x
		{
			for (int j = in2 - 1; j >= 3; j--) // y
			{
				if (i == j)
					res += grid[in1*j + i];
				else
					if (grid[in1*j + i] > 0)
						for (int k = 1; k <= j && k <= i - j; k++)
							grid[in1*k + i - j] += grid[in1*j + i];
			}

			if (grid[in1 * 1 + i] > 0)
				res += grid[in1 * 1 + i];
			if (grid[in1 * 2 + i] > 0)
				res += (i / 2)*grid[in1 * 2 + i];
		}

		OutSt.push((long double)res);

		delete[] grid;
	}
	else if (ofnum == 53) /// TODO: Floating point precision error, change to integer calc
	{
		long double n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		long double val = 0, tmp;
		bool isNeg = false;

		for (int i = 0; i <= n2; i++)
		{
			tmp = pow(n2 - i, n1);
			tmp /= tgamma(i + 1);
			tmp /= tgamma(n2 - i + 1);

			if (isNeg)
				val -= tmp;
			else
				val += tmp;

			isNeg = !isNeg;
		}

		OutSt.push(round(val));
	}
	else if (ofnum == 54)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(tgamma(n));
	}
	/*else if (ofnum == 55) /// TODO: E(n)
	{
		long double n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push();
	}*/

	return;
}

int OPrecComp(char o1, char o2)
{
	if ((o1 == 4 && o2 == 5) || (o1 == 5 && o2 == 4) || (o1 == 6 && o2 == 7) || (o1 == 7 && o2 == 6))
		return 0;
	else if (o1 < o2) // precedence: o1 > o2
		return 1;
	else
		return -1;
}

char FHash(char str[])
{
	if (strcmp(str, "max") == 0)
		return 21;
	else if (strcmp(str, "min") == 0)
		return 22;

	else if (strcmp(str, "sin") == 0)
		return 23;
	else if (strcmp(str, "cos") == 0)
		return 24;
	else if (strcmp(str, "tan") == 0)
		return 25;

	else if (strcmp(str, "csc") == 0)
		return 26;
	else if (strcmp(str, "sec") == 0)
		return 27;
	else if (strcmp(str, "cot") == 0)
		return 28;

	else if (strcmp(str, "asin") == 0)
		return 29;
	else if (strcmp(str, "acos") == 0)
		return 30;
	else if (strcmp(str, "atan") == 0)
		return 31;

	else if (strcmp(str, "acsc") == 0)
		return 32;
	else if (strcmp(str, "asec") == 0)
		return 33;
	else if (strcmp(str, "acot") == 0)
		return 34;

	else if (strcmp(str, "sinh") == 0)
		return 35;
	else if (strcmp(str, "cosh") == 0)
		return 36;
	else if (strcmp(str, "tanh") == 0)
		return 37;

	else if (strcmp(str, "csch") == 0)
		return 38;
	else if (strcmp(str, "sech") == 0)
		return 39;
	else if (strcmp(str, "coth") == 0)
		return 40;

	else if (strcmp(str, "asinh") == 0)
		return 41;
	else if (strcmp(str, "acosh") == 0)
		return 42;
	else if (strcmp(str, "atanh") == 0)
		return 43;

	else if (strcmp(str, "acsch") == 0)
		return 44;
	else if (strcmp(str, "asech") == 0)
		return 45;
	else if (strcmp(str, "acoth") == 0)
		return 46;

	else if (strcmp(str, "log") == 0)
		return 47;

	else if (strcmp(str, "erf") == 0)
		return 48;
	else if (strcmp(str, "erfc") == 0)
		return 49;

	else if (strcmp(str, "Perm") == 0)
		return 50;
	else if (strcmp(str, "Comb") == 0)
		return 51;
	else if (strcmp(str, "P") == 0)
		return 52;
	else if (strcmp(str, "S") == 0)
		return 53;

	else if (strcmp(str, "gamma") == 0)
		return 54;

	else if (strcmp(str, "E") == 0)
		return 55;

	return -1;
}

// with functions with variable argument, need to separate somehow with indicators

/*
Operator/Function Hash Values:
1 = LParen
2 = Fact / unary
3 = Exp / R
4 = Mult / L
5 = Div / L
6 = Add / L
7 = Sub / L
8~20 = Reserved
(1>2>3>(4, 5)>(6, 7))

=========================================

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

47 = log  //(base, exponand)

48 = erf
49 = erfc

50 = Perm
51 = Comb
52 = P (Partition of integer #)
53 = S (Subset partition)

54 = gamma

55 = E (Complete elliptic integral of the second kind)

~255 reserved
*/

/*
Mathematical constants are used with a preceding '#'
Lists: #pi, #e, #phi
*/