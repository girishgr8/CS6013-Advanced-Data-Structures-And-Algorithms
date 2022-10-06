
/* 											CS6013: ADVANCED DATA STRUCTURES AND ALGORITHMS (ADSA)
											Programming Assignment I - MULTIPLYING TWO POLYNOMIALS (Using Naive Approach and other using FFT)

NAME : 			THATTE GIRISH MAKARAND
ROLL NUMBER : 	CS22MTECH11005
 
*/

import java.util.*;

class main_prog {
	public static void main(String[] args) {
		Scanner sc = new Scanner(System.in);
		int[] poly_1 = new int[16];
		int[] poly_2 = new int[16];
		int[] naive_prod = new int[32];
		int[] fft_prod = new int[32];
		int deg_1, deg_2;

		// Taking input the degree of 1st polynomial into deg_1
		System.out.printf("Enter the degree of the first polynomial: ");
		deg_1 = sc.nextInt();
		sc.nextLine();
		// Taking input the coefficients of 1st polynomial
		System.out.printf(
				"Enter the " + (deg_1 + 1)
						+ " coefficients of the first polynomial in the increasing order of the degree of the monomials they belong to: ");
		for (int i = 0; i <= deg_1; i++)
			poly_1[i] = sc.nextInt();
		sc.nextLine();
		// Taking input the degree of 2nd polynomial into deg_2
		System.out.printf("Enter the degree of the second polynomial: ");
		deg_2 = sc.nextInt();
		sc.nextLine();
		// Taking input the coefficients of 2nd polynomial
		System.out.printf("Enter the " + (deg_2 + 1)
				+ " coefficients of the second polynomial in the increasing order of the degree of the monomials they belong to: ");
		for (int i = 0; i <= deg_2; i++)
			poly_2[i] = sc.nextInt();

		// Printing the 1st polynomial in human-readable format
		System.out.println("The first polynomial is: ");
		print_poly(poly_1, deg_1);

		// Printing the 2nd polynomial in human-readable format
		System.out.println("The second polynomial is: ");
		print_poly(poly_2, deg_2);

		// Performing the polynomial multiplication using naive approach in O(n^2)
		System.out.println("The product of the two polynomials obtained via naive polynomial multiplication is:");
		naive_prod = naive_polynomial_multiplication(poly_1, poly_2);

		// Printing the polynomial multiplication performed using naive approach
		print_poly(naive_prod, deg_1 + deg_2 + 1);

		// find the smallest power of 2 greater than deg_1 + deg_2 + 1 as N
		int N = find_N(deg_1, deg_2);

		// Performing the polynomial multiplication using FFT approach in O(n*logn)
		// Evaluating the 1st polynomial at N points i.e. in the point-value
		// representation form
		double[][] poly_1_evaluations = Eval(poly_1, N);
		// Evaluating the 2nd polynomial at N points i.e. in the point-value
		// representation form
		double[][] poly_2_evaluations = Eval(poly_2, N);

		// Multiplying the evaluated polynomials
		double[][] product_poly_evaluations = product_polynomial_evaluations(poly_1_evaluations, poly_2_evaluations, N);

		// Interpolating the point-value representation to coefficient form
		// representation
		double[][] result = EvalIFFT(product_poly_evaluations, N);

		System.out.println("The product of the two polynomials obtained via polynomial multiplication using FFT is:");
		for (int i = 0; i < N; i++)
			fft_prod[i] = (int) Math.round(result[i][0]);

		// Printing out the resultant polynomial generated using FFT
		print_poly(fft_prod, N);
		sc.close();
	}

	// naive_polynomial_multiplication() function : It does naive polynomial
	// multiplication in O(n^2) time
	public static int[] naive_polynomial_multiplication(int[] poly_1, int[] poly_2) {
		int[] naive_prod = new int[32];
		for (int i = 0; i < poly_1.length; i++) {
			for (int j = 0; j < poly_2.length; j++) {
				naive_prod[i + j] += poly_1[i] * poly_2[j];
			}
		}
		return naive_prod;
	}

	// add_complex_num() function takes four parameters which are real and imaginary
	// parts of two complex numbers
	// 1st complex number : a + ib
	// 2nd complex number : c + id
	// returns the real and imaginary part of the resultant complex number
	public static double[] add_complex_num(double a, double b, double c, double d) {
		double[] sum = new double[2];
		sum[0] = a + c;
		sum[1] = b + d;
		return sum;
	}

	// multiply_complex_num() function takes four parameters which are real and
	// imaginary parts of two complex numbers
	// 1st complex number : a + ib
	// 2nd complex number : c + id
	// returns the real and imaginary part of the resultant complex number
	public static double[] multiply_complex_num(double a, double b, double c, double d) {
		double[] mult = new double[2];
		mult[0] = (a * c) - (b * d);
		mult[1] = (a * d) + (b * c);
		return mult;
	}

	// find_N() function takes two parameters deg_1 and deg_2 and returns the
	// smallest power of 2 greater than equal to (deg_1 + deg_2 + 1)
	public static int find_N(int deg_1, int deg_2) {
		int N = 0, n = deg_1 + deg_2 + 1;
		while (true) {
			if ((n / (1 << N)) > 0)
				N++;
			else
				break;
		}
		return (1 << N);
	}

	// find_complex_roots() finds the Nth roots of unity as follows:
	// [w^0, w^1, w^2, w^3, w^4, .................., w^(n -1)]
	public static double[][] find_complex_roots(int N) {
		double[][] roots = new double[N][2];
		for (int i = 0; i < N; i++) {
			double angle_in_radians = (2 * Math.PI * i) / (double) N;
			roots[i][0] = Math.cos(angle_in_radians);
			roots[i][1] = Math.sin(angle_in_radians);
		}
		return roots;
	}

	// Eval() function : It converts the polynomial from coefficient form to
	// point-value form
	public static double[][] Eval(int[] poly, int N) {
		double[][] poly_eval = new double[N][2];

		// find the Nth complex roots of unity which will be used further in evaluation
		// of polynomial
		double[][] roots = find_complex_roots(N);

		// Base case: when N = 1, evaluation at single point i.e. single coefficient
		if (N == 1) {
			poly_eval[0][0] = poly[0];
			poly_eval[0][1] = 0;
			return poly_eval;
		}

		int[] poly_even = new int[N / 2];
		int[] poly_odd = new int[N / 2];

		// Divide the polynomial's coefficient into even and odd term coefficients using
		// "Divide and Conquer" strategy
		for (int i = 0; i < N / 2; i++) {
			poly_even[i] = poly[2 * i];
			poly_odd[i] = poly[2 * i + 1];
		}

		// Recursively evaluate for even term coefficients of polynomial
		double[][] poly_eval_even = Eval(poly_even, N / 2);
		// Recursively evaluate for odd term coefficients of polynomial
		double[][] poly_eval_odd = Eval(poly_odd, N / 2);

		for (int i = 0; i < N / 2; i++) {

			// poly_eval[i] = poly_eval_even[i] + (wn)^i * poly_eval_odd[i]
			double[] mult = multiply_complex_num(roots[i][0], roots[i][1], poly_eval_odd[i][0],
					poly_eval_odd[i][1]);

			double[] add_left = add_complex_num(poly_eval_even[i][0], poly_eval_even[i][1], mult[0], mult[1]);
			poly_eval[i][0] = add_left[0];
			poly_eval[i][1] = add_left[1];

			double[] add_right = add_complex_num(poly_eval_even[i][0], poly_eval_even[i][1], (-mult[0]), (-mult[1]));

			// poly_eval[i] = poly_eval_even[i] - (wn)^i * poly_eval_odd[i]
			poly_eval[i + N / 2][0] = add_right[0];
			poly_eval[i + N / 2][1] = add_right[1];
		}
		// return the evaluated polynomial at N points
		return poly_eval;
	}

	// Multiplying the two polynomials A(x) and B(x) in the point-value
	// representation as : C(x) = A(x) * B(x)
	public static double[][] product_polynomial_evaluations(double[][] poly_1_evaluations,
			double[][] poly_2_evaluations, int N) {
		double[][] product_poly_evaluations = new double[N][2];
		for (int i = 0; i < N; i++) {
			double a = poly_1_evaluations[i][0];
			double b = poly_1_evaluations[i][1];
			double c = poly_2_evaluations[i][0];
			double d = poly_2_evaluations[i][1];
			double[] mult = multiply_complex_num(a, b, c, d);

			product_poly_evaluations[i][0] = mult[0];
			product_poly_evaluations[i][1] = mult[1];
		}
		return product_poly_evaluations;
	}

	// EvalIFFT() function : It converts the polynomial from point-value form
	// coefficient form
	public static double[][] EvalIFFT(double[][] product_poly_evaluations, int N) {
		double[][] poly_eval = new double[N][2];

		// find the Nth complex roots of unity which will be used further in evaluation
		// of polynomial
		double[][] roots = find_complex_roots(N);

		// Base case: when N = 1, evaluation at single point i.e. single coefficient
		if (N == 1) {
			poly_eval[0][0] = product_poly_evaluations[0][0];
			poly_eval[0][1] = product_poly_evaluations[0][1];
			return poly_eval;
		}

		double[][] poly_even = new double[N / 2][2];
		double[][] poly_odd = new double[N / 2][2];

		// Divide the polynomial's coefficient into even and odd term coefficients using
		// "Divide and Conquer" strategy
		for (int i = 0; i < N / 2; i++) {
			poly_even[i][0] = product_poly_evaluations[2 * i][0];
			poly_even[i][1] = product_poly_evaluations[2 * i][1];

			poly_odd[i][0] = product_poly_evaluations[2 * i + 1][0];
			poly_odd[i][1] = product_poly_evaluations[2 * i + 1][1];
		}

		// Recursively evaluate for even term coefficients of polynomial
		double[][] poly_eval_even = EvalIFFT(poly_even, N / 2);
		// Recursively evaluate for odd term coefficients of polynomial
		double[][] poly_eval_odd = EvalIFFT(poly_odd, N / 2);

		for (int i = 0; i < N / 2; i++) {

			// poly_eval[i] = poly_eval_even[i] + (wn)^i * poly_eval_odd[i]
			double[] mult = multiply_complex_num(roots[i][0], -roots[i][1], poly_eval_odd[i][0],
					poly_eval_odd[i][1]);

			double[] add_left = add_complex_num(poly_eval_even[i][0], poly_eval_even[i][1], mult[0], mult[1]);
			poly_eval[i][0] = add_left[0] / 2;
			poly_eval[i][1] = add_left[1] / 2;

			// poly_eval[i] = poly_eval_even[i] - (wn)^i * poly_eval_odd[i]
			double[] add_right = add_complex_num(poly_eval_even[i][0], poly_eval_even[i][1], (-mult[0]), (-mult[1]));
			poly_eval[i + N / 2][0] = add_right[0] / 2;
			poly_eval[i + N / 2][1] = add_right[1] / 2;
		}
		// return the evaluated polynomial at N points
		return poly_eval;
	}

	// print_poly() function prints the polynomial in human-readable and easy format
	public static void print_poly(int[] poly, int deg) {
		String result = "";
		for (int i = deg; i >= 0; i--) {
			if (i != 0 && poly[i] != 0) {
				if (poly[i] != 1)
					result += "" + Math.abs(poly[i]);
				result += "x*" + i;
			}
			if (i > 0 && poly[i - 1] != 0)
				result += (poly[i - 1] > 0) ? " + " : " - ";
			if (i == 0) {
				if (poly[i] != 0)
					result += "" + Math.abs(poly[i]);
			}
		}
		result = result.trim();
		if (result.charAt(0) == '+')
			result = result.substring(1);
		if (result.charAt(result.length() - 1) == '+' || result.charAt(result.length() - 1) == '-')
			result = result.substring(0, result.length() - 1);
		result = result.trim();
		System.out.println(result);
	}
}