////**********************************************************************
////  Author: Bill Hallahan
////  Date: January 30, 2003
////**********************************************************************

//using System;

//namespace PolynomialRootFinder
//{
//    class Polynomial
//    {
//        double[] coefficients;
//        int degree;

//        public Polynomial()
//        {
//            coefficients = null;
//        }

//        public Polynomial(double scalar)
//        {
//            coefficients = new double[1];
//            degree = 0;
//            coefficients[0] = scalar;
//        }

//        public Polynomial(double x_coefficient, double scalar)
//        {
//            coefficients = new double[2];
//            degree = 1;
//            SetToFirstOrderPolynomial(x_coefficient, scalar);
//        }


//        public Polynomial(double x_squared_coefficient, double x_coefficient, double scalar)
//        {
//            coefficients = new double[3];
//            degree = 2;
//            SetToQuadraticPolynomial(x_squared_coefficient, x_coefficient, scalar);
//        }

//        public Polynomial(double[] coefficients)
//        {
//            this.coefficients = new double[coefficients.Length];
//            degree = coefficients.Length - 1;
//            SetCoefficients(coefficients);
//        }

//        public Polynomial(Polynomial polynomial) : this(polynomial.coefficients)
//        {
//        }

//        public Polynomial(int degree)
//        {
//            coefficients = new double[degree + 1];
//            this.degree = degree;
//        }

//        void SetToFirstOrderPolynomial(double x_coefficient, double scalar)
//        {
//            double[] coefficient_array = new double[2];
//            coefficient_array[0] = scalar;
//            coefficient_array[1] = x_coefficient;
//            SetCoefficients(coefficient_array);
//        }

//        void SetToQuadraticPolynomial(double x_squared_coefficient, double x_coefficient, double scalar)
//        {
//            double[] coefficient_array = new double[3];
//            coefficient_array[0] = scalar;
//            coefficient_array[1] = x_coefficient;
//            coefficient_array[2] = x_squared_coefficient;
//            SetCoefficients(coefficient_array);
//        }

//        void SetCoefficients(double[] coefficient)
//        {
//            for (int i = 0; i < coefficient.Length; ++i)
//                this.coefficients[i] = coefficient[i];
//        }


//        //======================================================================
//        //  Member Function: EvaluateReal
//        //
//        //  Abstract:
//        //
//        //    This method evaluates the polynomial using the passed real
//        //    value x. The algorithm used is Horner's method.
//        //
//        //
//        //  Input:
//        //
//        //    xr    A real value.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns a value of type 'double' that is equal
//        //    to the polynomial evaluated at the passed value x. 
//        //
//        //======================================================================

//        double EvaluateReal(double xr)
//        {
//            double pr = coefficients[degree];
//            int i = 0;

//            for (i = degree - 1; i >= 0; --i)
//            {
//                pr = pr * xr + coefficients[i];
//            }

//            return pr;
//        }

//        //======================================================================
//        //  Member Function: EvaluateReal
//        //
//        //  Abstract:
//        //
//        //    This method evaluates the polynomial using the passed real
//        //    value x. The algorithm used is Horner's method.
//        //
//        //
//        //  Input:
//        //
//        //    xr    A real value.
//        //
//        //    dr    A reference to a double which contains the real term
//        //          that the polynomial derivative evaluates to.
//        //
//        //  Return Value:
//        //
//        //    This function returns a value of type 'double' that is equal
//        //    to the polynomial evaluated at the passed value x. 
//        //
//        //======================================================================

//        double EvaluateReal(double xr, out double dr)
//        {
//            double pr = coefficients[degree];
//            dr = pr;

//            int i = 0;

//            for (i = degree - 1; i > 0; --i)
//            {
//                pr = pr * xr + coefficients[i];
//                dr = dr * xr + pr;
//            }

//            pr = pr * xr + coefficients[0];

//            return pr;
//        }

//        //======================================================================
//        //  Member Function: EvaluateImaginary
//        //
//        //  Abstract:
//        //
//        //    This function evaluates the value of a polynomial for a
//        //    specified pure imaginary value xi. The polynomial value
//        //    is evaluated by Horner's method.
//        //
//        //
//        //  Input:
//        //
//        //    xi        A double which equals the imaginary term used to
//        //              evaluate the polynomial.
//        //
//        //  Outputs:
//        //
//        //    pr        A reference to a double which contains the real term
//        //              that the polynomial evaluates to.
//        //
//        //    pi        A reference to a double which contains the imaginary
//        //              term that the polynomial evaluates to.
//        //
//        //  Return Value:
//        //
//        //    This function has no return value. 
//        //
//        //======================================================================

//        void EvaluateImaginary(double xi, out double pr, out double pi)
//        {
//            pr = coefficients[degree];
//            pi = 0;

//            int i = 0;

//            for (i = degree - 1; i >= 0; --i)
//            {
//                double temp = -pi * xi + coefficients[i];
//                pi = pr * xi;
//                pr = temp;
//            }

//            return;
//        }

//        //======================================================================
//        //  Member Function: EvaluateComplex
//        //
//        //  Abstract:
//        //
//        //    This function evaluates the value of a polynomial for a
//        //    specified complex value xr + i xi. The polynomial value
//        //    is evaluated by Horner's method.
//        //
//        //
//        //  Input:
//        //
//        //    xr        A double which equals the real term used to evaluate
//        //              the polynomial.
//        //
//        //    xi        A double which equals the imaginary term used to
//        //              evaluate the polynomial.
//        //
//        //  Outputs:
//        //
//        //    pr        A reference to a double which contains the real term
//        //              that the polynomial evaluates to.
//        //
//        //    pi        A reference to a double which contains the imaginary
//        //              term that the polynomial evaluates to.
//        //
//        //  Return Value:
//        //
//        //    This function has no return value. 
//        //
//        //======================================================================

//        void EvaluateComplex(double xr, double xi, out double pr, double pi)
//        {
//            pr = coefficients[degree];
//            pi = 0;

//            int i = 0;

//            for (i = degree - 1; i >= 0; --i)
//            {
//                double temp = pr * xr - pi * xi + coefficients[i];
//                pi = pr * xi + pi * xr;
//                pr = temp;
//            }

//            return;
//        }

//        //======================================================================
//        //  Member Function: EvaluateComplex
//        //
//        //  Abstract:
//        //
//        //    This function evaluates the value of a polynomial and the
//        //    value of the polynomials derivative for a specified complex
//        //    value xr + i xi. The polynomial value is evaluated by
//        //    Horner's method. The combination of the polynomial evaluation
//        //    and the derivative evaluation is known as the Birge-Vieta method.
//        //
//        //
//        //  Input:
//        //
//        //    xr        A double which equals the real term used to evaluate
//        //              the polynomial.
//        //
//        //    xi        A double which equals the imaginary term used to
//        //              evaluate the polynomial.
//        //
//        //  Outputs:
//        //
//        //    pr        A reference to a double which contains the real term
//        //              that the polynomial evaluates to.
//        //
//        //    pi        A reference to a double which contains the imaginary
//        //              term that the polynomial evaluates to.
//        //
//        //    dr        A reference to a double which contains the real term
//        //              that the polynomial derivative evaluates to.
//        //
//        //    di        A reference to a double which contains the imaginary
//        //              term that the polynomial derivative evaluates to.
//        //
//        //  Return Value:
//        //
//        //    This function has no return value. 
//        //
//        //======================================================================

//        void EvaluateComplex(double xr, double xi, out double pr, out double pi, out double dr, out double di)
//        {
//            pr = coefficients[degree];
//            pi = 0;
//            dr = pr;
//            di = 0;

//            double temp = 0.0;
//            int i = 0;

//            for (i = degree - 1; i > 0; --i)
//            {
//                temp = pr * xr - pi * xi + coefficients[i];
//                pi = pr * xi + pi * xr;
//                pr = temp;

//                temp = dr * xr - di * xi + pr;
//                di = dr * xi + di * xr + pi;
//                dr = temp;
//            }

//            temp = pr * xr - pi * xi + coefficients[0];
//            pi = pr * xi + pi * xr;
//            pr = temp;

//            return;
//        }

//        //======================================================================
//        //  Member Function: Derivative
//        //
//        //  Abstract:
//        //
//        //    This method calculates the derivative polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    None
//        //
//        //  Return Value:
//        //
//        //    This function returns a polynomial that is the derivative
//        //    of this polynomial.
//        //
//        //======================================================================

//        Polynomial Derivative()
//        {
//            Polynomial derivative_polynomial = new Polynomial(degree - 1);

//            //------------------------------------------------------------------
//            //  If this polynomial is just a scalar, i.e. it is of degree
//            //  zero then the derivative is zero.
//            //------------------------------------------------------------------

//            if (degree > 0)
//            {
//                //--------------------------------------------------------------
//                //  Calculate the derivative polynomial.
//                //--------------------------------------------------------------

//                double[] temp_ptr = derivative_polynomial.coefficients;

//                for (int i = degree; i > 0; --i)
//                {
//                    temp_ptr[i - 1] = (double)(i) * coefficients[i];
//                }
//            }

//            return derivative_polynomial;
//        }

//        //======================================================================
//        //  Member Function: Integral
//        //
//        //  Abstract:
//        //
//        //    This method calculates the integral polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    None
//        //
//        //  Return Value:
//        //
//        //    This function returns a polynomial that is the integral
//        //    of this polynomial.
//        //
//        //======================================================================

//        Polynomial Integral()
//        {
//            Polynomial integral_polynomial = new Polynomial(degree + 1);

//            //------------------------------------------------------------------
//            //  Calculate the integral polynomial.
//            //------------------------------------------------------------------

//            double[] temp_ptr = integral_polynomial.coefficients;
//            int i = 0;

//            for (i = degree; i > 0; --i)
//            {
//                temp_ptr[i + 1] = coefficients[i] / (double)(i + 1);
//            }

//            return integral_polynomial;
//        }

//        //======================================================================
//        //  Member Function: Degree
//        //
//        //  Abstract:
//        //
//        //    This method gets the polynomial degree.
//        //
//        //
//        //  Input:
//        //
//        //    None.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns a value of type 'int' that is the
//        //    degree of the polynomial.
//        //
//        //======================================================================

//        int Degree()
//        {
//            return degree;
//        }

//        //======================================================================
//        //  Member Function: FindRoots
//        //
//        //  Abstract:
//        //
//        //    This method determines the roots of a polynomial which has
//        //    real coefficients.
//        //
//        //
//        //  Input:
//        //
//        //
//        //    real_zero_vector_ptr       A double precision vector that will
//        //                               contain the real parts of the roots
//        //                               of the polynomial when this method
//        //                               returns.
//        //
//        //    imaginary_zero_vector_ptr  A double precision vector that will
//        //                               contain the real parts of the roots
//        //                               of the polynomial when this method
//        //                               returns.
//        //
//        //    roots_found_ptr           A pointer to an integer that will
//        //                              equal the number of roots found when
//        //                              this method returns. If the method
//        //                              returns SUCCESS then this value will
//        //                              always equal the degree of the
//        //                              polynomial.
//        //
//        //  Return Value:
//        //
//        //    This function returns an enum value of type
//        //    'RootStatus_T'.
//        //
//        //======================================================================

//        PolynomialRootFinder.RootStatus_T FindRoots(double[] real_zero_vector_ptr, double[] imaginary_zero_vector_ptr, int[] roots_found_ptr)
//        {

//            PolynomialRootFinder  polynomial_root_finder_ptr = new PolynomialRootFinder();

//            PolynomialRootFinder RootStatus_T status = root_finder_ptr.FindRoots(coefficients, real_zero_vector_ptr, imaginary_zero_vector_ptr, roots_found_ptr);
//            return status;
//        }

//        //======================================================================
//        //  Member Function: IncludeRealRoot
//        //
//        //  Abstract:
//        //
//        //    This method multiplies this polynomial by a first order term
//        //    that incorporates the real root.
//        //
//        //
//        //  Input:
//        //
//        //    real_value    A real root value.
//        //
//        //
//        //  Return Value:
//        //
//        //    The function has no return value.
//        //
//        //======================================================================

//        void IncludeRealRoot(double real_value)
//        {
//            double[] coefficient_array =new double[2];
//            coefficient_array[0] = -real_value;
//            coefficient_array[1] = 1.0;
//            Polynomial temp_polynomial = new Polynomial(coefficient_array);
//            this *= temp_polynomial;
//            return;
//        }

//        //======================================================================
//        //  Member Function: IncludeComplexConjugateRootPair
//        //
//        //  Abstract:
//        //
//        //    This method multiplies this polynomial by a second order
//        //    polynomial that incorporates a complex root pair.
//        //
//        //
//        //  Input:
//        //
//        //    real_value    A real root value.
//        //
//        //    imag_value    An imaginary root value.
//        //
//        //
//        //  Return Value:
//        //
//        //    The function has no return value.
//        //
//        //======================================================================

//        void IncludeComplexConjugateRootPair(double real_value, double imag_value)
//        {
//            double[] coefficient_array= new double[3];
//            coefficient_array[0] = real_value * real_value + imag_value * imag_value;
//            coefficient_array[1] = -(real_value + real_value);
//            coefficient_array[2] = 1.0;
//            Polynomial temp_polynomial = new Polynomial(coefficient_array);
//            this *= temp_polynomial;
//        }

//        //======================================================================
//        //  Member Function: Divide
//        //
//        //  Abstract:
//        //
//        //    This method divides this polynomial by a passed divisor
//        //    polynomial. The remainder polynomial is stored in the passed
//        //    instance remainder_polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    divisor_polynomial      The divisor polynomial
//        //
//        //    quotient_polynomial     A reference to an instance of class
//        //                            Polynomial that will contain the quotient
//        //                            polynomial when this method returns.
//        //
//        //    remainder_polynomial    A reference to an instance of class
//        //                            Polynomial that will contain the remainder
//        //                            polynomial when this method returns.
//        //
//        //  Return Value:
//        //
//        //    This function returns a value of type 'bool' that false if this
//        //    method fails. This can only occur if the divisor polynomial is
//        //    equal to the scalar value zero. Otherwise the return value is
//        //    true.
//        //
//        //======================================================================

//        bool Divide(Polynomial divisor_polynomial, out Polynomial quotient_polynomial, out Polynomial remainder_polynomial)
//        {
//            //------------------------------------------------------------------
//            //  If the divisor is zero then fail.
//            //------------------------------------------------------------------

//            int divisor_degree = divisor_polynomial.Degree();

//            bool non_zero_divisor_flag = (divisor_polynomial.Degree() >= 0);

//            if (non_zero_divisor_flag)
//            {
//                //--------------------------------------------------------------
//                //  If this dividend polynomial's degree is not greater than
//                //  or equal to the divisor polynomial's degree then do the division.
//                //--------------------------------------------------------------

//                remainder_polynomial = new Polynomial(this);
//                int dividend_degree = Degree();
//                quotient_polynomial = null;
//                int quotient_maximum_degree = dividend_degree - divisor_degree + 1;
//                quotient_polynomial = new Polynomial(quotient_maximum_degree);
//                quotient_polynomial.degree = -1;
//                double[] quotient_coefficient_ptr = quotient_polynomial.coefficients;
//                double[] dividend_coefficient_ptr = remainder_polynomial.coefficients;
//                double leading_divisor_coefficient = divisor_polynomial[divisor_degree];

//                //--------------------------------------------------------------
//                //  Loop and subtract each scaled divisor polynomial
//                //  to perform the division.
//                //--------------------------------------------------------------

//                int dividend_index = dividend_degree;

//                for (dividend_index = dividend_degree;
//                     dividend_index >= divisor_degree;
//                     --dividend_index)
//                {
//                    //----------------------------------------------------------
//                    //  Subtract the scaled divisor from the dividend.
//                    //----------------------------------------------------------

//                    double scale_value = remainder_polynomial[dividend_index] / leading_divisor_coefficient;

//                    //----------------------------------------------------------
//                    //  Increase the order of the quotient polynomial.
//                    //----------------------------------------------------------

//                    quotient_polynomial.degree += 1;
//                    int j = 0;

//                    for (j = quotient_polynomial.degree; j >= 1; --j)
//                    {
//                        quotient_coefficient_ptr[j] = quotient_coefficient_ptr[j - 1];
//                    }

//                    quotient_coefficient_ptr[0] = scale_value;

//                    //----------------------------------------------------------
//                    //  Subtract the scaled divisor from the dividend.
//                    //----------------------------------------------------------

//                    int dividend_degree_index = dividend_index;

//                    for (j = divisor_degree; j >= 0; --j)
//                    {
//                        dividend_coefficient_ptr[dividend_degree_index] -= divisor_polynomial[j] * scale_value;
//                        --dividend_degree_index;
//                    }
//                }

//                //--------------------------------------------------------------
//                //  Adjust the order of the current dividend polynomial.
//                //  This is the remainder polynomial.
//                //--------------------------------------------------------------

//                remainder_polynomial.AdjustPolynomialDegree();
//                quotient_polynomial.AdjustPolynomialDegree();
//            }
//            else
//            {
//                quotient_polynomial = new Polynomial(-1);
//                remainder_polynomial = new Polynomial(-1);
//            }

//            return non_zero_divisor_flag;
//        }
        
//        void AdjustPolynomialDegree()
//        {
//            while (degree > 0 && (Math.Abs(coefficients[degree]) < 2.2204460492503131e-016))
//            {
//                coefficients[degree] = 0.0;
//                degree--;
//            }

//            return;
//        }

//        //======================================================================
//        //  Member Function: operator []
//        //
//        //  Abstract:
//        //
//        //    This method returns the specified polynomial coefficient.
//        //
//        //
//        //  Input:
//        //
//        //    power_index      The coefficient index.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns a value of type 'double' that is the
//        //    coefficient value corresponding to the specified power.
//        //
//        //======================================================================

//        double this[int power_index]
//        {
//            //------------------------------------------------------------------
//            //  Ensure that the index is within range.
//            //------------------------------------------------------------------
//            get
//            {
//                if ((power_index < 0) || (power_index > degree))
//                {
//                    throw new Exception("Polynomial index out of range");
//                }

//                return coefficients[power_index];
//            }
//        }

//        //======================================================================
//        //  Member Function: operator +=
//        //
//        //  Abstract:
//        //
//        //    This method adds a polynomial to this polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    polynomial    An instance of class Polynomial
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        static public Polynomial operator +(Polynomial polynomial1, Polynomial polynomial2)
//        {
//            int i = 0;

//            if (polynomial1.degree >= polynomial2.degree)
//            {
//                for (i = 0; i <= polynomial2.degree; ++i)
//                {
//                    polynomial1.coefficients[i] += polynomial2.coefficients[i];
//                }
//            }
//            else
//            {

//                for (i = 0; i <= degree; ++i)
//                {
//                    polynomial1.coefficients[i] += polynomial2.coefficients[i];
//                }

//                for (i = degree + 1; i <= polynomial1.degree; ++i)
//                {
//                    polynomial1.coefficients[i] = polynomial2.coefficients[i];
//                }

//                degree = polynomial1.degree;
//            }

//            //------------------------------------------------------------------
//            //  If the leading coefficient(s) are zero, then decrease the
//            //  polynomial degree.
//            //------------------------------------------------------------------

//            AdjustPolynomialDegree();

//            return *this;
//        }

//        //======================================================================
//        //  Member Function: operator +=
//        //
//        //  Abstract:
//        //
//        //    This method adds a scalar to this polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    scalar    A scalar value.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        static public Polynomial operator +(Polynomial polynomial, double scalar)
//        {
//            assert(m_degree >= 0);

//            polynomial.coefficients[0] += scalar;

//            return this;
//        }

//        //======================================================================
//        //  Member Function: operator -=
//        //
//        //  Abstract:
//        //
//        //    This method subtracts a polynomial from this polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    polynomial    An instance of class Polynomial
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        static public Polynomial operator -(Polynomial polynomial)
//        {
//            assert(m_degree >= 0);

//            int i = 0;

//            if (m_degree >= polynomial.m_degree)
//            {
//                for (i = 0; i <= polynomial.m_degree; ++i)
//                {
//                    coefficients[i] -= polynomial.coefficients[i];
//                }
//            }
//            else
//            {
//                SetLength(polynomial.m_degree + 1, true);

//                for (i = 0; i <= m_degree; ++i)
//                {
//                    coefficients[i] -= polynomial.coefficients[i];
//                }

//                for (i = m_degree + 1; i <= polynomial.m_degree; ++i)
//                {
//                    coefficients[i] = -polynomial.coefficients[i];
//                }

//                m_degree = polynomial.m_degree;
//            }

//            //------------------------------------------------------------------
//            //  If the leading coefficient(s) are zero, then decrease the
//            //  polynomial degree.
//            //------------------------------------------------------------------

//            AdjustPolynomialDegree();

//            return *this;
//        }

//        //======================================================================
//        //  Member Function: operator -=
//        //
//        //  Abstract:
//        //
//        //    This method subtracts a scalar from this polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    scalar    A scalar value.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        static public Polynomial operator -(Polynomial polynomial, double scalar)
//        {
//            assert(m_degree >= 0);

//            polynomial.coefficients[0] -= scalar;

//            return this;
//        }

//        //======================================================================
//        //  Member Function: operator *=
//        //
//        //  Abstract:
//        //
//        //    This method multiplies a polynomial times this polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    polynomial    An instance of class Polynomial
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        static public operator *(Polynomial polynomial1, Polynomial polynomial2)
//{
//    //------------------------------------------------------------------
//    //  Create a temporary buffer to hold the product of the two
//    //  polynomials.
//    //------------------------------------------------------------------

//    assert(m_degree >= 0);

//    int convolution_length = m_degree + polynomial.m_degree + 1;

//    std::vector<double> temp_vector;
//    temp_vector.resize(convolution_length + 1);
//    double * temp_vector_ptr = &temp_vector[0];

//    //------------------------------------------------------------------
//    //  Zero the temporary buffer.
//    //------------------------------------------------------------------

//    int i = 0;

//    for (i = 0; i < convolution_length; ++i)
//    {
//        temp_vector_ptr[i] = 0.0;
//    }

//    //------------------------------------------------------------------
//    //  Calculate the convolution in the temporary buffer.
//    //------------------------------------------------------------------

//    for (i = 0; i <= m_degree; ++i)
//    {
//        for (int j = 0; j <= polynomial1.m_degree; ++j)
//        {
//            temp_vector_ptr[i + j] += polynomial1.coefficients[i] * polynomial2.coefficients[j];
//        }
//    }

//    //------------------------------------------------------------------
//    //  Make sure this buffer is large enough for the product.
//    //------------------------------------------------------------------

//    SetLength((unsigned int)(convolution_length), false);

//    //------------------------------------------------------------------
//    //  Store the result in this instance.
//    //------------------------------------------------------------------

//    m_degree = convolution_length - 1;

//    for (i = 0; i <= m_degree; ++i)
//    {
//        coefficients[i] = temp_vector_ptr[i];
//    }

//    //------------------------------------------------------------------
//    //  If the leading coefficient(s) are zero, then decrease the
//    //  polynomial degree.
//    //------------------------------------------------------------------

//    AdjustPolynomialDegree();

//    return *this;
//}

//        //======================================================================
//        //  Member Function: operator *=
//        //
//        //  Abstract:
//        //
//        //    This method multiplies a scalar time this polynomial.
//        //
//        //
//        //  Input:
//        //
//        //    scalar    A scalar value.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        static public Polynomial operator *(Polynomial polynomial, double scalar)
//        {
//            assert(m_degree >= 0);

//            int i = 0;

//            for (i = 0; i <= m_degree; ++i)
//            {
//                polynomal.m_coefficient_vector_ptr[i] *= scalar;
//            }

//            //------------------------------------------------------------------
//            //  If the leading coefficient(s) are zero, then decrease the
//            //  polynomial degree.
//            //------------------------------------------------------------------

//            AdjustPolynomialDegree();

//            return *this;
//        }

//        //======================================================================
//        //  Member Function: operator /=
//        //
//        //  Abstract:
//        //
//        //    This method divides this polynomial by a scalar.
//        //
//        //
//        //  Input:
//        //
//        //    scalar    A scalar value.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        static public Polynomial operator /(Polynomial polynomial, double scalar)
//        {
//            assert(m_degree >= 0);

//            int i = 0;

//            for (i = 0; i <= m_degree; ++i)
//            {
//                pollynomial.m_coefficient_vector_ptr[i] /= scalar;
//            }

//            return *this;
//        }

//        //======================================================================
//        //  Member Function: operator +
//        //
//        //  Abstract:
//        //
//        //    This method implements unary operator +()
//        //
//        //
//        //  Input:
//        //
//        //    None.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns a polynomial which is equal to this instance.
//        //
//        //======================================================================

//        object Clone()
//        {
//            assert(m_degree >= 0);
//            return this;
//        }

//        //======================================================================
//        //  Member Function: operator -
//        //
//        //  Abstract:
//        //
//        //    This method implements unary operator -().
//        //    This polynomials coefficients became the negative of
//        //    their present value and then this polynomial is returned.
//        //
//        //
//        //  Input:
//        //
//        //    None.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns a polynomial which is the negative of
//        //    this instance.
//        //
//        //======================================================================


//        //======================================================================
//        //  Member Function: operator =
//        //
//        //  Abstract:
//        //
//        //    This method sets this polynomial to a scalar value.
//        //
//        //
//        //  Input:
//        //
//        //    scalar    A scalar value.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        Polynomial set(double scalar)
//        {
//            SetCoefficients(scalar, 0);
//            return this;
//        }

//        //======================================================================
//        //  Member Function: AdjustPolynomialDegree
//        //
//        //  Abstract:
//        //
//        //    This method will decrease the polynomial degree until leading
//        //    coefficient is non-zero or until the polynomial degree is zero.
//        //
//        //
//        //  Input:
//        //
//        //    None.
//        //
//        //
//        //  Return Value:
//        //
//        //    This method has no return value.
//        //
//        //======================================================================

//        //======================================================================
//        //  Member Function: Copy
//        //
//        //  Abstract:
//        //
//        //    This method copies a passed polynomial into this instance.
//        //
//        //
//        //  Input:
//        //
//        //    polynomial    An instance of class Polynomial.
//        //
//        //
//        //  Return Value:
//        //
//        //    This function returns this polynomial.
//        //
//        //======================================================================

//        void Copy(Polynomial polynomial)
//        {
//            for (int i = 0; i <= polynomial.coefficients.Length; ++i)
//            {
//                coefficients[i] = polynomial.coefficients[i];
//            }

//            return;
//        }

//        //======================================================================
//        //  Global operators
//        //======================================================================

//        //======================================================================
//        //  Addition of an instance of the Polynomial class and a scalar.
//        //======================================================================

//        static public Polynomial operator +(double scalar, Polynomial polynomial)
//        {
//            return Polynomial(polynomial) += scalar;
//        }

//        //======================================================================
//        //  Subtraction of two instances of this class.
//        //======================================================================

//        static public Polynomial operator -(Polynomial minuend_polynomial, Polynomial subtrahend_polynomial)
//        {
//            return Polynomial(minuend_polynomial) -= subtrahend_polynomial;
//        }

//        //======================================================================
//        //  Subtraction with an instance of the Polynomial class and a scalar.
//        //======================================================================

//        static public Polynomial operator -(double scalar, Polynomial polynomial)
//        {
//            return (-Polynomial(polynomial)) + scalar;
//        }

//        //======================================================================
//        //  Multiplication of two instances of this class.
//        //======================================================================


//        //======================================================================
//        //  Multiplication of an instance of the Polynomial class and a scalar.
//        //======================================================================


//        //======================================================================
//        //  Division with an instance of the Polynomial class and a scalar.
//        //======================================================================

//    }
//}

