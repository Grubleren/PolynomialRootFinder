    ////**********************************************************************
////  Author: Bill Hallahan
////  Date: January 30, 2003
////**********************************************************************

using JH.Calculations;
using System;

namespace PolynomialRootFinder
{
    public class PolynomialRootFinder
    {
        public enum RootStatus
        {
            SUCCESS,
            LEADING_COEFFICIENT_IS_ZERO,
            SCALAR_VALUE_HAS_NO_ROOTS,
            FAILED_TO_CONVERGE
        };


        const float f_BASE = 2.0f;
        const float f_ETA = 1.192092896e-07F;
        const float f_ETA_N = (float)(10.0) * f_ETA;
        const float f_ETA_N_SQUARED = (float)(100.0) * f_ETA;
        const float f_MAXIMUM_FLOAT = float.MaxValue;
        const float f_MINIMUM_FLOAT = 1.175494351e-38F;
        const float f_XX_INITIAL_VALUE = (float)(0.70710678);
        const float f_COSR_INITIAL_VALUE = (float)(-0.069756474);
        const float f_SINR_INITIAL_VALUE = (float)(0.99756405);

        int m_degree;

        double[] m_p_vector_ptr;
        double[] m_qp_vector_ptr;
        double[] m_k_vector_ptr;
        double[] m_qk_vector_ptr;
        double[] m_svk_vector_ptr;

        int m_n;
        int m_n_plus_one;
        double m_real_s;
        double m_imag_s;
        double m_u;
        double m_v;
        double m_a = 0;
        double m_b = 0;
        double m_c = 0;
        double m_d = 0;
        double m_a1;
        double m_a3;
        double m_a7;
        double m_e;
        double m_f;
        double m_g;
        double m_h;
        double m_real_sz;
        double m_imag_sz;
        double m_real_lz = 0;
        double m_imag_lz = 0;
        double m_are;
        double m_mre;

        //======================================================================
        //  Member Function: FindRoots
        //
        //  Abstract:
        //
        //    This method determines the roots of a polynomial which
        //    has real coefficients. This code is based on FORTRAN
        //    code published in reference [1]. The method is based on
        //    an algorithm the three-stage algorithm described in
        //    Jenkins and Traub [2].
        //
        // 1. "Collected Algorithms from ACM, Volume III", Algorithms 493-545
        //    1983. (The root finding algorithms is number 493)
        //
        // 2. Jenkins, M. A. and Traub, J. F., "A three-stage algorithm for
        //    real polynomials using quadratic iteration", SIAM Journal of
        //    Numerical Analysis, 7 (1970), 545-566
        //
        // 3. Jenkins, M. A. and Traub, J. F., "Principles for testing
        //    polynomial zerofinding programs", ACM TOMS 1,
        //    1 (March 1975), 26-34
        //
        //
        //  Input:
        //
        //    All vectors below must be at least a length equal to degree + 1.
        //
        //    coefficicent_ptr           A double precision vector that contains
        //                               the polynomial coefficients in order
        //                               of increasing power.
        //
        //    degree                     The degree of the polynomial.
        //
        //    real_zero_vector_ptr       A double precision vector that will
        //                               contain the real parts of the roots
        //                               of the polynomial when this method
        //                               returns.
        //
        //    imaginary_zero_vector_ptr  A double precision vector that will
        //                               contain the real parts of the roots
        //                               of the polynomial when this method
        //                               returns.
        //
        //    number_of_roots_found_ptr  A pointer to an integer that will
        //                               equal the number of roots found when
        //                               this method returns. If the method
        //                               returns SUCCESS then this value will
        //                               always equal the degree of the
        //                               polynomial.
        //
        //  Return Value:
        //
        //    The function returns an enum value of type
        //    'RootStatus_T'.
        //
        //======================================================================

        public RootStatus FindRootsReverse(double[] coefficient_vector_ptr,
                               Complex[] zero_vector_ptr,
                               out int number_of_roots_found_ptr)
        {
            double[] coeff = new double[coefficient_vector_ptr.Length];
            double[] real_zero_vector_ptr = new double[zero_vector_ptr.Length];
            double[] imaginary_zero_vector_ptr = new double[zero_vector_ptr.Length];
            
            for (int i = 0; i < coeff.Length; i++)
                coeff[i] = coefficient_vector_ptr[coefficient_vector_ptr.Length - 1 - i];
            
            RootStatus status = FindRoots(coeff, coeff.Length - 1, real_zero_vector_ptr, imaginary_zero_vector_ptr, out number_of_roots_found_ptr);

            for (int i = 0; i < number_of_roots_found_ptr; i++)
            {
                zero_vector_ptr[i] = new Complex(real_zero_vector_ptr[i], imaginary_zero_vector_ptr[i]);
            }

            return status;
        }

        public RootStatus FindRoots(double[] coefficient_vector_ptr,
                               int degree,
                               double[] real_zero_vector_ptr,
                               double[] imaginary_zero_vector_ptr,
                               out int number_of_roots_found_ptr)
        {
            //--------------------------------------------------------------
            //  The algorithm fails if the polynomial is not at least
            //  degree on or the leading coefficient is zero.
            //--------------------------------------------------------------

            RootStatus status;
            number_of_roots_found_ptr = -1;
            if (degree == 0)
            {
                status = RootStatus.SCALAR_VALUE_HAS_NO_ROOTS;
            }
            else if (coefficient_vector_ptr[degree] == 0.0)
            {
                status = RootStatus.LEADING_COEFFICIENT_IS_ZERO;
            }
            else
            {
                //--------------------------------------------------------------
                //  Allocate temporary vectors used to find the roots.
                //--------------------------------------------------------------

                m_degree = degree;

                double[] temp_vector_ptr = new double[m_degree + 1];
                double[] pt_vector_ptr = new double[m_degree + 1];

                m_p_vector_ptr = new double[m_degree + 1];
                m_qp_vector_ptr = new double[m_degree + 1];
                m_k_vector_ptr = new double[m_degree + 1];
                m_qk_vector_ptr = new double[m_degree + 1];
                m_svk_vector_ptr = new double[m_degree + 1];

                //--------------------------------------------------------------
                //  m_are and m_mre refer to the unit error in + and *
                //  respectively. they are assumed to be the same as
                //  f_ETA.
                //--------------------------------------------------------------

                m_are = f_ETA;
                m_mre = f_ETA;
                double lo = f_MINIMUM_FLOAT / f_ETA;

                //--------------------------------------------------------------
                // Initialization of constants for shift rotation.
                //--------------------------------------------------------------

                double xx = f_XX_INITIAL_VALUE;
                double yy = -xx;
                double cosr = f_COSR_INITIAL_VALUE;
                double sinr = f_SINR_INITIAL_VALUE;
                m_n = m_degree;
                m_n_plus_one = m_n + 1;

                //--------------------------------------------------------------
                //  Make a copy of the coefficients in reverse order.
                //--------------------------------------------------------------

                int ii = 0;

                for (ii = 0; ii < m_n_plus_one; ++ii)
                {
                    m_p_vector_ptr[m_n - ii] = coefficient_vector_ptr[ii];
                }

                //--------------------------------------------------------------
                //  Assume failure. The status is set to SUCCESS if all
                //  the roots are found.
                //--------------------------------------------------------------

                status = RootStatus.FAILED_TO_CONVERGE;

                //--------------------------------------------------------------
                //  If there are any zeros at the origin, remove them.
                //--------------------------------------------------------------

                int jvar = 0;

                while (m_p_vector_ptr[m_n] == 0.0)
                {
                    jvar = m_degree - m_n;
                    real_zero_vector_ptr[jvar] = 0.0;
                    imaginary_zero_vector_ptr[jvar] = 0.0;
                    m_n_plus_one = m_n_plus_one - 1;
                    m_n = m_n - 1;
                }

                //--------------------------------------------------------------
                //  Loop and find polynomial zeros. In the original algorithm
                //  this loop was an infinite loop. Testing revealed that the
                //  number of main loop iterations to solve a polynomial of a
                //  particular degree is usually about half the degree.
                //  We loop twice that to make sure the solution is found.
                //  (This should be revisited as it might preclude solving
                //  some large polynomials.)
                //--------------------------------------------------------------

                int count = 0;

                for (count = 0; count < m_degree; ++count)
                {
                    //----------------------------------------------------------
                    //  Check for less than 2 zeros to finish the solution.
                    //----------------------------------------------------------

                    if (m_n <= 2)
                    {
                        if (m_n > 0)
                        {
                            //--------------------------------------------------
                            //  Calculate the final zero or pair of zeros.
                            //--------------------------------------------------

                            if (m_n == 1)
                            {
                                real_zero_vector_ptr[m_degree - 1] =
                                    -m_p_vector_ptr[1] / m_p_vector_ptr[0];

                                imaginary_zero_vector_ptr[m_degree - 1] = 0.0;
                            }
                            else
                            {
                                double x1;
                                double x2;
                                double x3;
                                double x4;

                                SolveQuadraticEquation(
                                    m_p_vector_ptr[0],
                                    m_p_vector_ptr[1],
                                    m_p_vector_ptr[2],
                                    out x1,
                                    out x2,
                                    out x3,
                                    out x4);

                                    real_zero_vector_ptr[m_degree - 2] = x1;
                                    imaginary_zero_vector_ptr[m_degree - 2]=x2;
                                    real_zero_vector_ptr[m_degree - 1]=x3;
                                    imaginary_zero_vector_ptr[m_degree - 1] = x4;
                            }
                        }

                        m_n = 0;
                        status = RootStatus.SUCCESS;
                        break;
                    }
                    else
                    {
                        //------------------------------------------------------
                        //  Find largest and smallest moduli of coefficients.
                        //------------------------------------------------------

                        double max = 0.0;
                        double min = f_MAXIMUM_FLOAT;
                        double xvar;

                        for (ii = 0; ii < m_n_plus_one; ++ii)
                        {
                            xvar = (double)(Math.Abs((double)(m_p_vector_ptr[ii])));

                            if (xvar > max)
                            {
                                max = xvar;
                            }

                            if ((xvar != 0.0) && (xvar < min))
                            {
                                min = xvar;
                            }
                        }

                        //------------------------------------------------------
                        //  Scale if there are large or very small coefficients.
                        //  Computes a scale factor to multiply the coefficients
                        //  of the polynomial. The scaling is done to avoid
                        //  overflow and to avoid undetected underflow from
                        //  interfering with the convergence criterion.
                        //  The factor is a power of the base.
                        //------------------------------------------------------

                        bool do_scaling_flag = false;
                        double sc = lo / min;

                        if (sc <= 1.0)
                        {
                            do_scaling_flag = f_MAXIMUM_FLOAT / sc < max;
                        }
                        else
                        {
                            do_scaling_flag = max < 10.0;

                            if (!do_scaling_flag)
                            {
                                if (sc == 0.0)
                                {
                                    sc = f_MINIMUM_FLOAT;
                                }
                            }
                        }

                        //------------------------------------------------------
                        //  Conditionally scale the data.
                        //------------------------------------------------------

                        if (do_scaling_flag)
                        {
                            int lvar = (int)(Math.Log(sc) / Math.Log(f_BASE) + 0.5);
                            double factor = Math.Pow((double)(f_BASE * 1.0), (double)lvar);

                            if (factor != 1.0)
                            {
                                for (ii = 0; ii < m_n_plus_one; ++ii)
                                {
                                    m_p_vector_ptr[ii] = factor * m_p_vector_ptr[ii];
                                }
                            }
                        }

                        //------------------------------------------------------
                        //  Compute lower bound on moduli of zeros.
                        //------------------------------------------------------

                        for (ii = 0; ii < m_n_plus_one; ++ii)
                        {
                            pt_vector_ptr[ii] = (double)(Math.Abs((double)(m_p_vector_ptr[ii])));
                        }

                        pt_vector_ptr[m_n] = -pt_vector_ptr[m_n];

                        //------------------------------------------------------
                        //  Compute upper estimate of bound.
                        //------------------------------------------------------

                        xvar = (double)
                            (Math.Exp((Math.Log(-pt_vector_ptr[m_n]) - Math.Log(pt_vector_ptr[0]))
                                / (double)(m_n)));

                        //------------------------------------------------------
                        //  If newton step at the origin is better, use it.
                        //------------------------------------------------------

                        double xm;

                        if (pt_vector_ptr[m_n - 1] != 0.0)
                        {
                            xm = -pt_vector_ptr[m_n] / pt_vector_ptr[m_n - 1];

                            if (xm < xvar)
                            {
                                xvar = xm;
                            }
                        }

                        //------------------------------------------------------
                        //  Chop the interval (0, xvar) until ff <= 0
                        //------------------------------------------------------

                        double ff;

                        while (true)
                        {
                            xm = (double)(xvar * 0.1);
                            ff = pt_vector_ptr[0];

                            for (ii = 1; ii < m_n_plus_one; ++ii)
                            {
                                ff = ff * xm + pt_vector_ptr[ii];
                            }

                            if (ff <= 0.0)
                            {
                                break;
                            }

                            xvar = xm;
                        }

                        double dx = xvar;

                        //------------------------------------------------------
                        //  Do newton iteration until xvar converges to two
                        //  decimal places.
                        //------------------------------------------------------

                        while (true)
                        {
                            if ((double)(Math.Abs(dx / xvar)) <= 0.005)
                            {
                                break;
                            }

                            ff = pt_vector_ptr[0];
                            double df = ff;

                            for (ii = 1; ii < m_n; ++ii)
                            {
                                ff = ff * xvar + pt_vector_ptr[ii];
                                df = df * xvar + ff;
                            }

                            ff = ff * xvar + pt_vector_ptr[m_n];
                            dx = ff / df;
                            xvar = xvar - dx;
                        }

                        double bnd = xvar;

                        //------------------------------------------------------
                        //  Compute the derivative as the intial m_k_vector_ptr
                        //  polynomial and do 5 steps with no shift.
                        //------------------------------------------------------

                        int n_minus_one = m_n - 1;

                        for (ii = 1; ii < m_n; ++ii)
                        {
                            m_k_vector_ptr[ii] =
                                (double)(m_n - ii) * m_p_vector_ptr[ii] / (double)(m_n);
                        }

                        m_k_vector_ptr[0] = m_p_vector_ptr[0];
                        double aa = m_p_vector_ptr[m_n];
                        double bb = m_p_vector_ptr[m_n - 1];
                        bool zerok_flag = m_k_vector_ptr[m_n - 1] == 0.0;

                        int jj = 0;

                        for (jj = 1; jj <= 5; ++jj)
                        {
                            double cc = m_k_vector_ptr[m_n - 1];

                            if (zerok_flag)
                            {
                                //----------------------------------------------
                                //  Use unscaled form of recurrence.
                                //----------------------------------------------

                                for (jvar = n_minus_one; jvar > 0; --jvar)
                                {
                                    m_k_vector_ptr[jvar] = m_k_vector_ptr[jvar - 1];
                                }

                                m_k_vector_ptr[0] = 0.0;
                                zerok_flag = m_k_vector_ptr[m_n - 1] == 0.0;
                            }
                            else
                            {
                                //----------------------------------------------
                                //  Use scaled form of recurrence if value
                                //  of m_k_vector_ptr at 0 is nonzero.
                                //----------------------------------------------

                                double tvar = -aa / cc;

                                for (jvar = n_minus_one; jvar > 0; --jvar)
                                {
                                    m_k_vector_ptr[jvar] =
                                        tvar * m_k_vector_ptr[jvar - 1] + m_p_vector_ptr[jvar];
                                }

                                m_k_vector_ptr[0] = m_p_vector_ptr[0];
                                zerok_flag =
                                    Math.Abs(m_k_vector_ptr[m_n - 1]) <= Math.Abs(bb) * f_ETA_N;
                            }
                        }

                        //------------------------------------------------------
                        //  Save m_k_vector_ptr for restarts with new shifts.
                        //------------------------------------------------------

                        for (ii = 0; ii < m_n; ++ii)
                        {
                            temp_vector_ptr[ii] = m_k_vector_ptr[ii];
                        }

                        //------------------------------------------------------
                        //  Loop to select the quadratic corresponding to
                        //   each new shift.
                        //------------------------------------------------------

                        int cnt = 0;

                        for (cnt = 1; cnt <= 20; ++cnt)
                        {
                            //--------------------------------------------------
                            //  Quadratic corresponds to a double shift to a
                            //  non-real point and its complex conjugate. The
                            //  point has modulus 'bnd' and amplitude rotated
                            //  by 94 degrees from the previous shift.
                            //--------------------------------------------------

                            double xxx = cosr * xx - sinr * yy;
                            yy = sinr * xx + cosr * yy;
                            xx = xxx;
                            m_real_s = bnd * xx;
                            m_imag_s = bnd * yy;
                            m_u = -2.0 * m_real_s;
                            m_v = bnd;

                            //--------------------------------------------------
                            //  Second stage calculation, fixed quadratic.
                            //  Variable nz will contain the number of
                            //   zeros found when function Fxshfr() returns.
                            //--------------------------------------------------

                            int nz = Fxshfr(20 * cnt);

                            if (nz != 0)
                            {
                                //----------------------------------------------
                                //  The second stage jumps directly to one of
                                //  the third stage iterations and returns here
                                //  if successful. Deflate the polynomial,
                                //  store the zero or zeros and return to the
                                //  main algorithm.
                                //----------------------------------------------

                                jvar = m_degree - m_n;
                                real_zero_vector_ptr[jvar] = m_real_sz;
                                imaginary_zero_vector_ptr[jvar] = m_imag_sz;
                                m_n_plus_one = m_n_plus_one - nz;
                                m_n = m_n_plus_one - 1;

                                for (ii = 0; ii < m_n_plus_one; ++ii)
                                {
                                    m_p_vector_ptr[ii] = m_qp_vector_ptr[ii];
                                }

                                if (nz != 1)
                                {
                                    real_zero_vector_ptr[jvar + 1] = m_real_lz;
                                    imaginary_zero_vector_ptr[jvar + 1] = m_imag_lz;
                                }

                                break;

                                //----------------------------------------------
                                //  If the iteration is unsuccessful another
                                //  quadratic is chosen after restoring
                                //  m_k_vector_ptr.
                                //----------------------------------------------
                            }

                            for (ii = 0; ii < m_n; ++ii)
                            {
                                m_k_vector_ptr[ii] = temp_vector_ptr[ii];
                            }
                        }
                    }
                }

                //--------------------------------------------------------------
                //  If no convergence with 20 shifts then adjust the degree
                //  for the number of roots found.
                //--------------------------------------------------------------

                number_of_roots_found_ptr = m_degree - m_n;
            }

            return status;
        }

        //======================================================================
        //  Computes up to l2var fixed shift m_k_vector_ptr polynomials,
        //  testing for convergence in the linear or quadratic
        //  case. initiates one of the variable shift
        //  iterations and returns with the number of zeros
        //  found.
        //
        //    l2var  An integer that is the limit of fixed shift steps.
        //
        //  Return Value:
        //    nz   An integer that is the number of zeros found.
        //======================================================================

        int Fxshfr(int l2var)
        {
            //------------------------------------------------------------------
            //  Evaluate polynomial by synthetic division.
            //------------------------------------------------------------------

            QuadraticSyntheticDivision(m_n_plus_one,
                                       m_u,
                                       m_v,
                                       m_p_vector_ptr,
                                       m_qp_vector_ptr,
                                       ref m_a,
                                       ref m_b);

            int itype = CalcSc();

            int nz = 0;
            float betav = 0.25f;
            float betas = 0.25f;
            float oss = (float)(m_real_s);
            float ovv = (float)(m_v);
            float ots = 0;
            float otv = 0;
            double ui = 0;
            double vi = 0;
            double svar;

            int jvar = 0;

            for (jvar = 1; jvar <= l2var; ++jvar)
            {
                //--------------------------------------------------------------
                //  Calculate next m_k_vector_ptr polynomial and estimate m_v.
                //--------------------------------------------------------------

                NextK(itype);
                itype = CalcSc();
                Newest(itype, ref ui, ref vi);
                float vv = (float)(vi);

                //--------------------------------------------------------------
                //  Estimate svar
                //--------------------------------------------------------------

                float ss = 0.0f;

                if (m_k_vector_ptr[m_n - 1] != 0.0)
                {
                    ss = (float)(-m_p_vector_ptr[m_n] / m_k_vector_ptr[m_n - 1]);
                }

                float tv = 1.0f;
                float ts = 1.0f;

                if ((jvar != 1) && (itype != 3))
                {
                    //----------------------------------------------------------
                    //  Compute relative measures of convergence of
                    //  svar and m_v sequences.
                    //----------------------------------------------------------

                    if (vv != 0.0)
                    {
                        tv = (float)(Math.Abs((vv - ovv) / vv));
                    }

                    if (ss != 0.0)
                    {
                        ts = (float)(Math.Abs((ss - oss) / ss));
                    }

                    //----------------------------------------------------------
                    //  If decreasing, multiply two most recent convergence
                    //  measures.
                    //----------------------------------------------------------

                    float tvv = 1.0f;

                    if (tv < otv)
                    {
                        tvv = tv * otv;
                    }

                    float tss = 1.0f;

                    if (ts < ots)
                    {
                        tss = ts * ots;
                    }

                    //----------------------------------------------------------
                    //  Compare with convergence criteria.
                    //----------------------------------------------------------

                    bool vpass_flag = tvv < betav;
                    bool spass_flag = tss < betas;

                    if (spass_flag || vpass_flag)
                    {
                        //------------------------------------------------------
                        //  At least one sequence has passed the convergence
                        //  test. Store variables before iterating.
                        //------------------------------------------------------

                        double svu = m_u;
                        double svv = m_v;
                        int ii = 0;

                        for (ii = 0; ii < m_n; ++ii)
                        {
                            m_svk_vector_ptr[ii] = m_k_vector_ptr[ii];
                        }

                        svar = ss;

                        //------------------------------------------------------
                        //  Choose iteration according to the fastest
                        //  converging sequence.
                        //------------------------------------------------------

                        bool vtry_flag = false;
                        bool stry_flag = false;
                        bool exit_outer_loop_flag = false;

                        bool start_with_real_iteration_flag =
                            (spass_flag && ((!vpass_flag) || (tss < tvv)));

                        do
                        {
                            if (!start_with_real_iteration_flag)
                            {
                                nz = QuadraticIteration(ui, vi);

                                if (nz > 0)
                                {
                                    exit_outer_loop_flag = true;
                                    break;
                                }

                                //----------------------------------------------
                                //  Quadratic iteration has failed. flag
                                //  that it has been tried and decrease
                                //  the convergence criterion.
                                //----------------------------------------------

                                vtry_flag = true;
                                betav = (float)(betav * 0.25);
                            }

                            //--------------------------------------------------
                            //  Try linear iteration if it has not been
                            //  tried and the svar sequence is converging.
                            //--------------------------------------------------

                            if (((!stry_flag) && spass_flag)
                                || start_with_real_iteration_flag)
                            {
                                if (!start_with_real_iteration_flag)
                                {
                                    for (ii = 0; ii < m_n; ++ii)
                                    {
                                        m_k_vector_ptr[ii] = m_svk_vector_ptr[ii];
                                    }
                                }
                                else
                                {
                                    start_with_real_iteration_flag = false;
                                }

                                int iflag = 0;

                                nz = RealIteration(ref svar, ref iflag);

                                if (nz > 0)
                                {
                                    exit_outer_loop_flag = true;
                                    break;
                                }

                                //----------------------------------------------
                                //  Linear iteration has failed. Flag that
                                //  it has been tried and decrease the
                                //  convergence criterion.
                                //----------------------------------------------

                                stry_flag = true;
                                betas = (float)(betas * 0.25);

                                if (iflag != 0)
                                {
                                    //------------------------------------------
                                    //  If linear iteration signals an almost
                                    //  double real zero attempt quadratic
                                    //  iteration.
                                    //------------------------------------------

                                    ui = -(svar + svar);
                                    vi = svar * svar;

                                    continue;
                                }
                            }

                            //--------------------------------------------------
                            //  Restore variables
                            //--------------------------------------------------

                            m_u = svu;
                            m_v = svv;

                            for (ii = 0; ii < m_n; ++ii)
                            {
                                m_k_vector_ptr[ii] = m_svk_vector_ptr[ii];
                            }

                            //----------------------------------------------
                            //  Try quadratic iteration if it has not been
                            //  tried and the m_v sequence is converging.
                            //----------------------------------------------
                        }
                        while (vpass_flag && (!vtry_flag));

                        if (exit_outer_loop_flag)
                        {
                            break;
                        }

                        //------------------------------------------------------
                        //  Recompute m_qp_vector_ptr and scalar values to
                        //  continue the second stage.
                        //------------------------------------------------------

                        QuadraticSyntheticDivision(m_n_plus_one,
                                                   m_u,
                                                   m_v,
                                                   m_p_vector_ptr,
                                                   m_qp_vector_ptr,
                                                   ref m_a,
                                                   ref m_b);

                        itype = CalcSc();
                    }
                }

                ovv = vv;
                oss = ss;
                otv = tv;
                ots = ts;
            }

            return nz;
        }

        //======================================================================
        //  Variable-shift m_k_vector_ptr-polynomial iteration for
        //  a quadratic factor converges only if the zeros are
        //  equimodular or nearly so.
        //
        //    uu  Coefficients of starting quadratic
        //    vv  Coefficients of starting quadratic
        //
        //  Return value:
        //    nz  The number of zeros found.
        //======================================================================

        int QuadraticIteration(double uu, double vv)
        {
            //------------------------------------------------------------------
            //  Main loop
            //------------------------------------------------------------------

            double ui = 0.0;
            double vi = 0.0;
            float omp = 0.0F;
            float relstp = 0.0F;
            int itype = 0;
            bool tried_flag = false;
            int jvar = 0;
            int nz = 0;
            m_u = uu;
            m_v = vv;

            while (true)
            {
                SolveQuadraticEquation(1.0,
                                       m_u,
                                       m_v,
                                       out m_real_sz,
                                       out m_imag_sz,
                                       out m_real_lz,
                                       out m_imag_lz);

                //--------------------------------------------------------------
                //  Return if roots of the quadratic are real and not close
                //  to multiple or nearly equal and  of opposite sign.
                //--------------------------------------------------------------

                if (Math.Abs(Math.Abs(m_real_sz) - Math.Abs(m_real_lz)) > 0.01 * Math.Abs(m_real_lz))
                {
                    break;
                }

                //--------------------------------------------------------------
                //  Evaluate polynomial by quadratic synthetic division.
                //------------------------------------------------------------------

                QuadraticSyntheticDivision(m_n_plus_one,
                                           m_u,
                                           m_v,
                                           m_p_vector_ptr,
                                           m_qp_vector_ptr,
                                           ref m_a,
                                           ref m_b);

                float mp = (float)(Math.Abs(m_a - m_real_sz * m_b) + Math.Abs(m_imag_sz * m_b));

                //--------------------------------------------------------------
                //  Compute a rigorous  bound on the rounding error in
                //  evaluting m_p_vector_ptr.
                //--------------------------------------------------------------

                float zm = (float)(Math.Sqrt((float)(Math.Abs((float)(m_v)))));
                float ee = (float)(2.0 * (float)(Math.Abs((float)(m_qp_vector_ptr[0]))));
                float tvar = (float)(-m_real_sz * m_b);
                int ii = 0;

                for (ii = 1; ii < m_n; ++ii)
                {
                    ee = ee * zm + (float)(Math.Abs((float)(m_qp_vector_ptr[ii])));
                }

                ee = ee * zm + (float)(Math.Abs((float)(m_a) + tvar));
                ee = (float)((5.0 * m_mre + 4.0 * m_are) * ee
                    - (5.0 * m_mre + 2.0 * m_are) * ((float)(Math.Abs((float)(m_a) + tvar)) + (float)(Math.Abs((float)(m_b))) * zm)
                        + 2.0 * m_are * (float)(Math.Abs(tvar)));

                //--------------------------------------------------------------
                //  Iteration has converged sufficiently if the polynomial
                //  value is less than 20 times this bound.
                //--------------------------------------------------------------

                if (mp <= 20.0 * ee)
                {
                    nz = 2;
                    break;
                }

                jvar = jvar + 1;

                //--------------------------------------------------------------
                //  Stop iteration after 20 steps.
                //--------------------------------------------------------------

                if (jvar > 20)
                {
                    break;
                }

                if ((jvar >= 2) && ((relstp <= 0.01)
                    && (mp >= omp) && (!tried_flag)))
                {
                    //----------------------------------------------------------
                    //  A cluster appears to be stalling the convergence.
                    //  Five fixed shift steps are taken with a m_u, m_v
                    //  close to the cluster.
                    //----------------------------------------------------------

                    if (relstp < f_ETA)
                    {
                        relstp = f_ETA;
                    }

                    relstp = (float)(Math.Sqrt(relstp));
                    m_u = m_u - m_u * relstp;
                    m_v = m_v + m_v * relstp;

                    QuadraticSyntheticDivision(m_n_plus_one,
                                               m_u,
                                               m_v,
                                               m_p_vector_ptr,
                                               m_qp_vector_ptr,
                                               ref m_a,
                                               ref m_b);

                    for (ii = 0; ii < 5; ++ii)
                    {
                        itype = CalcSc();
                        NextK(itype);
                    }

                    tried_flag = true;
                    jvar = 0;
                }

                omp = mp;

                //--------------------------------------------------------------
                //  Calculate next m_k_vector_ptr polynomial and
                //  new m_u and m_v.
                //--------------------------------------------------------------

                itype = CalcSc();
                NextK(itype);
                itype = CalcSc();
                Newest(itype, ref ui, ref vi);

                //--------------------------------------------------------------
                //  If vi is zero the iteration is not converging.
                //--------------------------------------------------------------

                if (vi == 0.0)
                {
                    break;
                }

                relstp = (float)(Math.Abs((vi - m_v) / vi));
                m_u = ui;
                m_v = vi;
            }

            return nz;
        }

        //======================================================================
        //  Variable-shift h polynomial iteration for a real zero.
        //
        //    sss      Starting iterate
        //    flag     Flag to indicate a pair of zeros near real axis.
        //
        //  Return Value:
        //     Number of zero found.
        //======================================================================

        int RealIteration(ref double sss, ref int flag)
        {
            //------------------------------------------------------------------
            //  Main loop
            //------------------------------------------------------------------

            double tvar = 0.0;
            float omp = 0.0F;
            int nz = 0;
            flag = 0;
            int jvar = 0;
            double svar = sss;

            while (true)
            {
                double pv = m_p_vector_ptr[0];

                //--------------------------------------------------------------
                //  Evaluate m_p_vector_ptr at svar
                //--------------------------------------------------------------

                m_qp_vector_ptr[0] = pv;
                int ii = 0;

                for (ii = 1; ii < m_n_plus_one; ++ii)
                {
                    pv = pv * svar + m_p_vector_ptr[ii];
                    m_qp_vector_ptr[ii] = pv;
                }

                float mp = (float)(Math.Abs(pv));

                //--------------------------------------------------------------
                //  Compute a rigorous bound on the error in evaluating p
                //--------------------------------------------------------------

                double ms = (double)(Math.Abs(svar));
                double ee = (m_mre / (m_are + m_mre)) * (double)(Math.Abs((double)(m_qp_vector_ptr[0])));

                for (ii = 1; ii < m_n_plus_one; ++ii)
                {
                    ee = ee * ms + (float)(Math.Abs((double)(m_qp_vector_ptr[ii])));
                }

                //--------------------------------------------------------------
                //  Iteration has converged sufficiently if the
                //  polynomial value is less than 20 times this bound.
                //--------------------------------------------------------------

                if (mp <= 20.0 * ((m_are + m_mre) * ee - m_mre * mp))
                {
                    nz = 1;
                    m_real_sz = svar;
                    m_imag_sz = 0.0;
                    break;
                }

                jvar = jvar + 1;

                //--------------------------------------------------------------
                //  Stop iteration after 10 steps.
                //--------------------------------------------------------------

                if (jvar > 10)
                {
                    break;
                }

                if ((jvar >= 2)
                    && ((Math.Abs(tvar) <= 0.001 * Math.Abs(svar - tvar))
                    && (mp > omp)))
                {
                    //----------------------------------------------------------
                    //  A cluster of zeros near the real axis has been
                    //  encountered. Return with flag set to initiate
                    //  a quadratic iteration.
                    //----------------------------------------------------------

                    flag = 1;
                    sss = svar;
                    break;
                }

                //--------------------------------------------------------------
                //  Return if the polynomial value has increased significantly.
                //--------------------------------------------------------------

                omp = mp;

                //--------------------------------------------------------------
                //  Compute t, the next polynomial, and the new iterate.
                //--------------------------------------------------------------

                double kv = m_k_vector_ptr[0];
                m_qk_vector_ptr[0] = kv;

                for (ii = 1; ii < m_n; ++ii)
                {
                    kv = kv * svar + m_k_vector_ptr[ii];
                    m_qk_vector_ptr[ii] = kv;
                }

                if (Math.Abs(kv) <= Math.Abs(m_k_vector_ptr[m_n - 1]) * f_ETA_N)
                {
                    m_k_vector_ptr[0] = 0.0;

                    for (ii = 1; ii < m_n; ++ii)
                    {
                        m_k_vector_ptr[ii] = m_qk_vector_ptr[ii - 1];
                    }
                }
                else
                {
                    //----------------------------------------------------------
                    //  Use the scaled form of the recurrence if the
                    //  value of m_k_vector_ptr at svar is non-zero.
                    //----------------------------------------------------------

                    tvar = -pv / kv;
                    m_k_vector_ptr[0] = m_qp_vector_ptr[0];

                    for (ii = 1; ii < m_n; ++ii)
                    {
                        m_k_vector_ptr[ii] = tvar * m_qk_vector_ptr[ii - 1] + m_qp_vector_ptr[ii];
                    }
                }

                //--------------------------------------------------------------
                //  Use unscaled form.
                //--------------------------------------------------------------

                kv = m_k_vector_ptr[0];

                for (ii = 1; ii < m_n; ++ii)
                {
                    kv = kv * svar + m_k_vector_ptr[ii];
                }

                tvar = 0.0;

                if (Math.Abs(kv) > Math.Abs(m_k_vector_ptr[m_n - 1]) * f_ETA_N)
                {
                    tvar = -pv / kv;
                }

                svar = svar + tvar;
            }

            return nz;
        }

        //======================================================================
        //  This routine calculates scalar quantities used to compute
        //  the next m_k_vector_ptr polynomial and new estimates of the
        //  quadratic coefficients.
        //
        //  Return Value:
        //    type  Integer variable set here indicating how the
        //    calculations are normalized to avoid overflow.
        //======================================================================

        int CalcSc()
        {
            //------------------------------------------------------------------
            //  Synthetic division of m_k_vector_ptr by the quadratic 1, m_u, m_v.
            //------------------------------------------------------------------

            QuadraticSyntheticDivision(m_n,
                                       m_u,
                                       m_v,
                                       m_k_vector_ptr,
                                       m_qk_vector_ptr,
                                       ref m_c,
                                       ref m_d);

            int itype = 0;

            if ((Math.Abs(m_c) <= Math.Abs(m_k_vector_ptr[m_n - 1]) * f_ETA_N_SQUARED)
                && (Math.Abs(m_d) <= Math.Abs(m_k_vector_ptr[m_n - 2]) * f_ETA_N_SQUARED))
            {
                //--------------------------------------------------------------
                //  itype == 3 Indicates the quadratic is almost a
                //  factor of m_k_vector_ptr.
                //--------------------------------------------------------------

                itype = 3;
            }
            else if (Math.Abs(m_d) >= Math.Abs(m_c))
            {
                //--------------------------------------------------------------
                //  itype == 2 Indicates that all formulas are divided by m_d.
                //--------------------------------------------------------------

                itype = 2;
                m_e = m_a / m_d;
                m_f = m_c / m_d;
                m_g = m_u * m_b;
                m_h = m_v * m_b;
                m_a3 = (m_a + m_g) * m_e + m_h * (m_b / m_d);
                m_a1 = m_b * m_f - m_a;
                m_a7 = (m_f + m_u) * m_a + m_h;
            }
            else
            {
                //--------------------------------------------------------------
                //  itype == 1 Indicates that all formulas are divided by m_c.
                //--------------------------------------------------------------

                itype = 1;
                m_e = m_a / m_c;
                m_f = m_d / m_c;
                m_g = m_u * m_e;
                m_h = m_v * m_b;
                m_a3 = m_a * m_e + (m_h / m_c + m_g) * m_b;
                m_a1 = m_b - m_a * (m_d / m_c);
                m_a7 = m_a + m_g * m_d + m_h * m_f;
            }

            return itype;
        }

        //======================================================================
        //  Computes the next k polynomials using scalars computed in CalcSc.
        //======================================================================

        void NextK(int itype)
        {
            int ii = 0;

            if (itype == 3)
            {
                //--------------------------------------------------------------
                //  Use unscaled form of the recurrence if type is 3.
                //--------------------------------------------------------------

                m_k_vector_ptr[0] = 0.0;
                m_k_vector_ptr[1] = 0.0;

                for (ii = 2; ii < m_n; ++ii)
                {
                    m_k_vector_ptr[ii] = m_qk_vector_ptr[ii - 2];
                }
            }
            else
            {
                double temp = m_a;

                if (itype == 1)
                {
                    temp = m_b;
                }

                if (Math.Abs(m_a1) <= Math.Abs(temp) * f_ETA_N)
                {
                    //----------------------------------------------------------
                    //  If m_a1 is nearly zero then use a special form of
                    //  the recurrence.
                    //----------------------------------------------------------

                    m_k_vector_ptr[0] = 0.0;
                    m_k_vector_ptr[1] = -m_a7 * m_qp_vector_ptr[0];

                    for (ii = 2; ii < m_n; ++ii)
                    {
                        m_k_vector_ptr[ii] = m_a3 * m_qk_vector_ptr[ii - 2] - m_a7 * m_qp_vector_ptr[ii - 1];
                    }
                }
                else
                {
                    //----------------------------------------------------------
                    //  Use scaled form of the recurrence.
                    //----------------------------------------------------------

                    m_a7 = m_a7 / m_a1;
                    m_a3 = m_a3 / m_a1;
                    m_k_vector_ptr[0] = m_qp_vector_ptr[0];
                    m_k_vector_ptr[1] = m_qp_vector_ptr[1] - m_a7 * m_qp_vector_ptr[0];

                    for (ii = 2; ii < m_n; ++ii)
                    {
                        m_k_vector_ptr[ii] =
                            m_a3 * m_qk_vector_ptr[ii - 2] - m_a7 * m_qp_vector_ptr[ii - 1] + m_qp_vector_ptr[ii];
                    }
                }
            }

            return;
        }

        //======================================================================
        //  Compute new estimates of the quadratic coefficients using the
        //  scalars computed in CalcSc.
        //======================================================================

        void Newest(int itype, ref double uu, ref double vv)
        {
            //------------------------------------------------------------------
            //  Use formulas appropriate to setting of itype.
            //------------------------------------------------------------------

            if (itype == 3)
            {
                //--------------------------------------------------------------
                //  If itype == 3 the quadratic is zeroed.
                //--------------------------------------------------------------

                uu = 0.0;
                vv = 0.0;
            }
            else
            {
                double a4;
                double a5;

                if (itype == 2)
                {
                    a4 = (m_a + m_g) * m_f + m_h;
                    a5 = (m_f + m_u) * m_c + m_v * m_d;
                }
                else
                {
                    a4 = m_a + m_u * m_b + m_h * m_f;
                    a5 = m_c + (m_u + m_v * m_f) * m_d;
                }

                //--------------------------------------------------------------
                //  Evaluate new quadratic coefficients.
                //--------------------------------------------------------------

                double b1 = -m_k_vector_ptr[m_n - 1] / m_p_vector_ptr[m_n];
                double b2 = -(m_k_vector_ptr[m_n - 2] + b1 * m_p_vector_ptr[m_n - 1]) / m_p_vector_ptr[m_n];
                double c1 = m_v * b2 * m_a1;
                double c2 = b1 * m_a7;
                double c3 = b1 * b1 * m_a3;
                double c4 = c1 - c2 - c3;
                double temp = a5 + b1 * a4 - c4;

                if (temp != 0.0)
                {
                    uu = m_u - (m_u * (c3 + c2) + m_v * (b1 * m_a1 + b2 * m_a7)) / temp;
                    vv = m_v * (1.0 + c4 / temp);
                }
            }

            return;
        }

        //======================================================================
        //  Divides p by the quadratic  1, u, v placing the quotient in q
        //  and the remainder in a,b
        //======================================================================

        void QuadraticSyntheticDivision(int n_plus_one,
                                                              double u,
                                                              double v,
                                                              double[] p_ptr,
                                                              double[] q_ptr,
                                                              ref double a,
                                                              ref double b)
        {
            b = p_ptr[0];
            q_ptr[0] = b;
            a = p_ptr[1] - u * b;
            q_ptr[1] = a;

            int ii = 0;

            for (ii = 2; ii < n_plus_one; ++ii)
            {
                double c = p_ptr[ii] - u * a - v * b;
                q_ptr[ii] = c;
                b = a;
                a = c;
            }

            return;
        }

        //======================================================================
        //                                          2
        //  Calculate the zeros of the quadratic a x + b x + c.
        //  the quadratic formula, modified to avoid overflow, is used to find
        //  the larger zero if the zeros are real and both zeros are complex.
        //  the smaller real zero is found directly from the product of the
        //  zeros c / a.
        //======================================================================

        void SolveQuadraticEquation(double a,
                                                          double b,
                                                          double c,
                                                          out double sr,
                                                          out double si,
                                                          out double lr,
                                                          out double li)
        {
            if (a == 0.0)
            {
                if (b != 0.0)
                {
                    sr = -c / b;
                }
                else
                {
                    sr = 0.0;
                }

                lr = 0.0;
                si = 0.0;
                li = 0.0;
            }
            else if (c == 0.0)
            {
                sr = 0.0;
                lr = -b / a;
                si = 0.0;
                li = 0.0;
            }
            else
            {
                //--------------------------------------------------------------
                //  Compute discriminant avoiding overflow.
                //--------------------------------------------------------------

                double d;
                double e;
                double bvar = b / 2.0;

                if (Math.Abs(bvar) < Math.Abs(c))
                {
                    if (c < 0.0)
                    {
                        e = -a;
                    }
                    else
                    {
                        e = a;
                    }

                    e = bvar * (bvar / Math.Abs(c)) - e;

                    d = Math.Sqrt(Math.Abs(e)) * Math.Sqrt(Math.Abs(c));
                }
                else
                {
                    e = 1.0 - (a / bvar) * (c / bvar);
                    d = Math.Sqrt(Math.Abs(e)) * Math.Abs(bvar);
                }

                if (e >= 0.0)
                {
                    //----------------------------------------------------------
                    //  Real zeros
                    //----------------------------------------------------------

                    if (bvar >= 0.0)
                    {
                        d = -d;
                    }

                    lr = (-bvar + d) / a;
                    sr = 0.0;

                    if (lr != 0.0)
                    {
                        sr = (c / lr) / a;
                    }

                    si = 0.0;
                    li = 0.0;
                }
                else
                {
                    //----------------------------------------------------------
                    //  Complex conjugate zeros
                    //----------------------------------------------------------

                    sr = -bvar / a;
                    lr = sr;
                    si = Math.Abs(d / a);
                    li = -si;
                }
            }

            return;
        }

        public double[] BuildPolynomial(Complex[] poles)
        {
            Complex[] polynomial = new Complex[poles.Length + 1];
            double[] coefficients = new double[poles.Length + 1];
            polynomial[0] = new Complex(1, 0);
            for (int i = 0; i < poles.Length; i++)
            {
                for (int j = polynomial.Length - 1; j > 0; j--)
                {
                    polynomial[j] -= polynomial[j - 1] * poles[i];
                }
            }
            for (int j = 0; j < polynomial.Length; j++)
                coefficients[j] = polynomial[j].re;

            return coefficients;

        }

        public double[] MulPolynomials(double[] p0, double[] p1)
        {
            double[] result = new double[p0.Length + p1.Length - 1];

            for (int j = 0; j < p0.Length; j++)
            {
                for (int i = 0; i < p1.Length; i++)
                    if (i + j < result.Length)
                        result[i + j] += p0[j] * p1[i];
            }

            return result;
        }
    }
}
