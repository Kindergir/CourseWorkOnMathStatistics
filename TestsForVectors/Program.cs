using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestsForVectors
{
    class Program
    {
        public static double[][] InverseMatrix(double[][] matrix)
        {
            if (matrix.Length != 0 && matrix.Length != matrix[0].Length)
                throw new ArgumentException("Non-square matrix");

            if (Math.Abs(CalculateMatrixDeterminant(matrix)) < 0.000001)
                throw new ArgumentException("Degenerate matrix");

            int N = matrix[0].Length;
            double[][] Rev = new double[N][];

            for (int i = 0; i < N; i++)
                Rev[i] = new double[N];

            for (int i = 0; i < N; i++)
                Rev[i][i] = 1.0;

            for (int k = 0; k < N; k++)
            {
                for (int i = 0; i < N; i++)
                {
                    if (i != k)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            if (j != k)
                                matrix[i][j] -= matrix[k][j] * matrix[i][k] / matrix[k][k];
                            Rev[i][j] -= Rev[k][j] * matrix[i][k] / matrix[k][k];
                        }
                        matrix[i][k] = 0.0;
                    }
                }

                for (int i = 0; i < N; i++)
                {
                    if (i != k)
                        matrix[k][i] /= matrix[k][k];
                    Rev[k][i] /= matrix[k][k];
                }
                matrix[k][k] = 1.0;
            }
            return Rev;
        }

        static void Main(string[] args)
        {
            double[] a = {1.0, 2, 3};
            double[][] b = { new []{ 10.0, 5.0 }, new[] { 4.0, 2.0 }};

            var c1 = InverseMatrix(b);

            for (int i = 0; i < c1.Length; ++i)
            {
                for (int j = 0; j < c1[i].Length; ++j)
                    Console.Write(c1[i][j] + " ");
                Console.WriteLine();
            }
        }
    }
}
