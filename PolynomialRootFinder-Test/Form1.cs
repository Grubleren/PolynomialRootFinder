using System;
using System.Windows.Forms;

namespace PolynomialRootFinder
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            double[] coefficients = new double[5];
            int degree = 3;
            double[] real_zero_vector_ptr = new double[4];
            double[] imaginary_zero_vector_ptr = new double[4];
            int number_of_roots_found_ptr;
            PolynomialRootFinder finder = new PolynomialRootFinder();

            // Highest degree , highest index

            double alpha = 0;
            double beta = 0.3;
            double gamma = 0.3;

            coefficients[3] = 1.0;
            coefficients[2] = -alpha;
            coefficients[1] = -(beta+gamma);
            coefficients[0] = gamma;

            PolynomialRootFinder.RootStatus status = finder.FindRoots(coefficients, degree, real_zero_vector_ptr, imaginary_zero_vector_ptr, out number_of_roots_found_ptr);

            for (int i = 0; i < degree; i++)
                Console.WriteLine("{0}    {1}", real_zero_vector_ptr[i], imaginary_zero_vector_ptr[i]);

                Close();

        }
    }
}
