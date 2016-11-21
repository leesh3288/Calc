#pragma warning(disable: 4146) // unsigned type unary minus operator
#pragma warning(disable: 4800) // int value force cast to bool
#pragma warning(disable: 4996) // deprecation

#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <stack>
#include <map>
#include <algorithm>
#include <mpirxx.h>
#include <mpreal.h>
#include <iostream>

using mpfr::mpreal;

#define MAXL_CONST 10
#define MAXL_FUNC 10


std::stack<mpreal> OutSt; // originally an output QUEUE, implemented as a stack for easy calculation (however, calculations must be done in reverse order)
std::stack<char> OpSt;
std::map<const char*, mpreal> CV; // ConstVal, saves constant values as desired precision

char equ[100000], outtype[10];
int decprec;
bool isSCI;

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
	printf("Input Precision (Significant Digits)\n >> ");
	scanf("%d", &decprec);
	printf("Output in Normalized Scientific Notation(SCI) or Decimal Form(DEC)\n >> ");
	scanf("%s", outtype);
	for (int i = 0; outtype[i] != '\0'; i++)
		outtype[i] = toupper(outtype[i]);
	if (strcmp(outtype, "SCI") == 0)
		isSCI = true;

	mpreal::set_default_prec(mpfr::digits2bits(decprec + 1));

	return;
}

void process()
{
	// const values init
	CV["pi"] = mpfr::const_pi();
	CV["e"] = mpfr::const_euler();
	CV["catalan"] = mpfr::const_catalan();
	CV["eps"] = std::numeric_limits<mpreal>::epsilon(); // the smallest value eps such that 1 + eps != eps (classic machine epsilon) | eps is NOT A CONSTANT

	int ite = 0, eqlen = (int)strlen(equ);

	while (ite < eqlen)
	{
		if (equ[ite] == ' ') // ignore spacing
			ite++;
		else if (equ[ite] >= '0' && equ[ite] <= '9') // number
		{
			int st = ite;
			mpreal inp;
			for (; equ[ite] >= '0'&&equ[ite] <= '9'; ite++);

			if (equ[ite] == '.')
				for (ite++; equ[ite] >= '0' && equ[ite] <= '9'; ite++);

			mpfr_strtofr(inp.mpfr_ptr(), equ+st, NULL, 10, mpreal::get_default_rnd());

			OutSt.push(inp);
		}
		else if ((equ[ite] >= 'a'&&equ[ite] <= 'z') || (equ[ite] >= 'A'&&equ[ite] <= 'Z') || equ[ite] == '_') // function & constants (input w/o #)
		{
			char fcs[MAXL_FUNC + 1];
			int cct;

			for (cct = 0; (equ[ite] >= '0'&&equ[ite] <= '9') || (equ[ite] >= 'a'&&equ[ite] <= 'z') || (equ[ite] >= 'A'&&equ[ite] <= 'Z') || equ[ite] == '_'; ite++, cct++)
				fcs[cct] = equ[ite];
			fcs[cct] = '\0';

			if (strcmp(fcs, "pi") == 0) // checks if it's a predefined constant
			{
				OutSt.push(CV["pi"]);
				continue;
			}
			else if (strcmp(fcs, "e") == 0)
			{
				OutSt.push(CV["e"]);
				continue;
			}
			else if (strcmp(fcs, "catalan") == 0)
			{
				OutSt.push(CV["catalan"]);
				continue;
			}

			OpSt.push(FHash(fcs));
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
				printf("Error(P): Left parenthesis not encountered, parsing failed.\n");
				return;
			}

			ite++;
		}
		else if (equ[ite] == '(') // left paren
		{
			OpSt.push(1);

			ite++;

			if (equ[ite] == '-') // solution for the negative number input -n: (-n) => (-1*n)
			{
				OutSt.push(mpreal(-1));

				OpSt.push(4); // multiplication operator '*'

				ite++;
			}
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
				printf("Error(P): Left parenthesis not encountered, parsing failed.\n");
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
				if (equ[ite + 1] == '!')
				{
					ite++;
					OFProc(8);
				}
				else
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
		if (OpSt.top() == 1) // if LParen is left over
		{
			printf("Error(P): Left parenthesis remaining (right parenthesis missing), parsing failed.\n");
			return;
		}
		OFProc(OpSt.top());
		OpSt.pop();
	}

	return;
}

void output()
{
	bool term = false;

	if (OutSt.size() != 1)
	{
		printf("Error(O): Size of output stack != 1\n");
		term = true;
	}
	if (!OpSt.empty())
	{
		printf("Error(O): Operator stack not empty.\n");
		term = true;
	}

	if (!term)
	{
		if (isSCI)
		{
			mpfr_out_str(stdout, 10, decprec, OutSt.top().mpfr_srcptr(), mpreal::get_default_rnd());
			printf("\n");
		}
		else
		{
			char form[15];
			sprintf(form, "%%.%dRNg\n", decprec);
			mpfr_printf(form, OutSt.top().mpfr_srcptr());
		}
	}

	return;
}

void OFProc(char ofnum) // operator/function processor
{
	// operator
	if (ofnum == 2)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		mpfr_fac_ui(n.mpfr_ptr(), n.toULong(mpreal::get_default_rnd()), mpreal::get_default_rnd());

		OutSt.push(n);
	}
	else if (ofnum == 3)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(pow(n1, n2));
	}
	else if (ofnum == 4)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1*n2);
	}
	else if (ofnum == 5)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1 / n2);
	}
	else if (ofnum == 6)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1 + n2);
	}
	else if (ofnum == 7)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(n1 - n2);
	}
	else if (ofnum == 8)
	{
		mpreal n;
		mpz_class tmp;
		n = OutSt.top();
		OutSt.pop();

		mpz_2fac_ui(tmp.get_mpz_t(), n.toULLong(mpreal::get_default_rnd()));

		OutSt.push(mpreal(tmp.get_mpz_t()));
	}

	// function
	else if (ofnum == 21)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(std::max(n1, n2));
	}
	else if (ofnum == 22)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(std::min(n1, n2));
	}
	else if (ofnum == 23)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(sin(n));
	}
	else if (ofnum == 24)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(cos(n));
	}
	else if (ofnum == 25)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(tan(n));
	}
	else if (ofnum == 26)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / sin(n));
	}
	else if (ofnum == 27)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / cos(n));
	}
	else if (ofnum == 28)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / tan(n));
	}
	else if (ofnum == 29)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asin(n));
	}
	else if (ofnum == 30)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acos(n));
	}
	else if (ofnum == 31)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atan(n));
	}
	else if (ofnum == 32)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asin(1 / n));
	}
	else if (ofnum == 33)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acos(1 / n));
	}
	else if (ofnum == 34)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atan(1 / n));
	}
	else if (ofnum == 35)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(sinh(n));
	}
	else if (ofnum == 36)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(cosh(n));
	}
	else if (ofnum == 37)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(tanh(n));
	}
	else if (ofnum == 38)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / sinh(n));
	}
	else if (ofnum == 39)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / cosh(n));
	}
	else if (ofnum == 40)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(1 / tanh(n));
	}
	else if (ofnum == 41)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asinh(n));
	}
	else if (ofnum == 42)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acosh(n));
	}
	else if (ofnum == 43)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atanh(n));
	}
	else if (ofnum == 44)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(asinh(1 / n));
	}
	else if (ofnum == 45)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(acosh(1 / n));
	}
	else if (ofnum == 46)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(atanh(1 / n));
	}
	else if (ofnum == 47)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		OutSt.push(log(n2) / log(n1));
	}
	else if (ofnum == 48)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(erf(n));
	}
	else if (ofnum == 49)
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(erfc(n));
	}
	else if (ofnum == 50)
	{
		mpreal n1, n2;
		mpz_class tmp, t1, t2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		mpfr_get_z(t1.get_mpz_t(), n1.mpfr_srcptr(), mpreal::get_default_rnd());
		mpfr_get_z(t2.get_mpz_t(), n2.mpfr_srcptr(), mpreal::get_default_rnd());
		mpz_bin_ui(tmp.get_mpz_t(), t1.get_mpz_t(), t2.get_ui());
		mpz_fac_ui(t1.get_mpz_t(), t2.get_ui());

		OutSt.push(mpreal(((mpz_class)(tmp*t1)).get_mpz_t()));
	}
	else if (ofnum == 51)
	{
		mpreal n1, n2;
		mpz_class tmp, t1, t2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		mpfr_get_z(t1.get_mpz_t(), n1.mpfr_srcptr(), mpreal::get_default_rnd());
		mpfr_get_z(t2.get_mpz_t(), n2.mpfr_srcptr(), mpreal::get_default_rnd());
		mpz_bin_ui(tmp.get_mpz_t(), t1.get_mpz_t(), t2.get_ui());

		OutSt.push(mpreal(tmp.get_mpz_t()));
	}
	else if (ofnum == 52)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		int in1 = n1.toLong(mpreal::get_default_rnd()), in2 = n2.toLong(mpreal::get_default_rnd()); // input restricted to (signed) int range due to indexing & memory
		mpz_class res = 0;

		in1++; in2++;
		mpz_class *grid = new mpz_class[in1*in1]; // access by grid[in1*y+x]

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

		OutSt.push(mpreal(res.get_mpz_t()));

		delete[] grid;
	}
	else if (ofnum == 53)
	{
		mpreal n1, n2;
		mpz_class i1, i2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();
		mpfr_get_z(i1.get_mpz_t(), n1.mpfr_srcptr(), mpreal::get_default_rnd());
		mpfr_get_z(i2.get_mpz_t(), n2.mpfr_srcptr(), mpreal::get_default_rnd());

		mpz_class ft;
		mpreal val = 0, tmp;
		bool isNeg = (bool)(i2.get_ui()&1);

		for (mpz_class i = 0; i <= i2; i++)
		{
			mpz_pow_ui(ft.get_mpz_t(), i.get_mpz_t(), i1.get_ui());
			tmp = mpreal(ft.get_mpz_t());

			mpz_fac_ui(ft.get_mpz_t(), i.get_ui()); // input restricted to unsigned (long) int range due to internal factorial function parameter type
			tmp /= ft.get_mpz_t();
			mpz_fac_ui(ft.get_mpz_t(), ((mpz_class)(i2 - i)).get_ui());
			tmp /= ft.get_mpz_t();

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
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		OutSt.push(tgamma(n));
	}
	else if (ofnum == 55) // EllipticK(n), calculated by AGM
	{
		mpreal n;
		n = OutSt.top();
		OutSt.pop();

		mpfr_agm(n.mpfr_ptr(), mpreal((1 + n)).mpfr_srcptr(), mpreal((1 - n)).mpfr_srcptr(), mpreal::get_default_rnd());

		OutSt.push(CV["pi"] / (2 * n));
	}
	else if (ofnum == 56) // EllipticE(n), calculated by AGM, Landen Transformation & much more calculations
	{
		mpreal n, res = 0;
		n = OutSt.top();
		OutSt.pop();

		mpreal a0 = 0, b0 = 0, a1 = 1, b1 = sqrt(1 - n*n), c, powtwo = 1;
		n.setPrecision(mpreal::get_default_prec() + 5); // add about 5 precision (to make c converge less than epsilon)
		res.setPrecision(mpreal::get_default_prec() + 5);
		a0.setPrecision(mpreal::get_default_prec() + 5);
		b0.setPrecision(mpreal::get_default_prec() + 5);
		a1.setPrecision(mpreal::get_default_prec() + 5);
		b1.setPrecision(mpreal::get_default_prec() + 5);
		c.setPrecision(mpreal::get_default_prec() + 5);
		powtwo.setPrecision(mpreal::get_default_prec() + 5);

		c = sqrt(a1*a1 - b1*b1); // only for first term

		res += powtwo*c*c;

		while (abs(c) > mpfr::machine_epsilon(mpreal::get_default_prec() + 2)) // AGM | add about 2 precision & take machine epsilon
		{
			a0 = a1;
			b0 = b1;

			a1 = (a0 + b0) / 2;
			b1 = sqrt(a0*b0);

			c = (a0 - b0) / 2;
			powtwo *= 2;

			res += powtwo*c*c;
		}

		mpreal ans = (1 - res / 2)*CV["pi"] / (2 * a1);
		ans.setPrecision(mpreal::get_default_prec());
		OutSt.push(ans);
	}
	else if (ofnum == 57)
	{
		mpreal n1, n2;
		n2 = OutSt.top();
		OutSt.pop();
		n1 = OutSt.top();
		OutSt.pop();

		mpfr_agm(n1.mpfr_ptr(), n1.mpfr_srcptr(), n2.mpfr_srcptr(), mpreal::get_default_rnd());

		OutSt.push(n1);
	}
	else if (ofnum == 58)
	{
		mpreal n;
		mpz_class tmp;
		n = OutSt.top();
		OutSt.pop();

		mpz_fib_ui(tmp.get_mpz_t(), n.toULLong(mpreal::get_default_rnd()));

		OutSt.push(mpreal(tmp.get_mpz_t()));
	}

	return;
}

int OPrecComp(char o1, char o2)
{
	if ((o1 == 4 && o2 == 5) || (o1 == 5 && o2 == 4) || (o1 == 6 && o2 == 7) || (o1 == 7 && o2 == 6) || (o1 == o2))
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

	else if (strcmp(str, "EllipticK") == 0)
		return 55;
	else if (strcmp(str, "EllipticE") == 0)
		return 56;

	else if (strcmp(str, "AGM") == 0)
		return 57;

	else if (strcmp(str, "fib") == 0)
		return 58;

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
8 = DoubleFact / unary
9~20 = Reserved
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

47 = log(base, exponand)

48 = erf
49 = erfc

50 = Perm
51 = Comb
52 = P (Partition of integer #)
53 = S (Subset partition)

54 = gamma

55 = EllipticK (Complete elliptic integral of the first kind)
56 = EllipticE (Complete elliptic integral of the second kind)

57 = AGM (Arithmetic-Geometric Mean)

58 = fib (nth fibonacci sequence)

~255 reserved
*/

/*
Constants list
Lists: pi, e, catalan
*/