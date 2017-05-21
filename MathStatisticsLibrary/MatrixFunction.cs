using System;
using System.Collections.Generic;

namespace MathStatisticsLibrary
{
    public class MatrixFunction
    {
        public static double[] TransposeVector(double[] vector)
        {
            double[] result = new double[vector.Length];
            for (int i = 0, j = vector.Length - 1; i < vector.Length; ++i, --j)
                result[i] = vector[j];
            return result;
        }

        public static double AlgebraicComplement(double[][] matrix, int row, int column)
        {
            double[][] subMatrix = DeleteRowAndColumnInMatrix(matrix, row, column);
            return Math.Pow(-1, (row + column + 2)) * MatrixDeterminant(subMatrix);
        }

        public static double MatrixDeterminant(double[][] matrix)
        {
            if (matrix.Length == 0 || matrix.Length != matrix[0].Length)
                throw new ArgumentException("For this matrix there is no determinant");

            double[][] l, u;
            DecomposeIntoLU(matrix, out l, out u);
            return UpperTriangularMatrixDeterminant(u);
        }

        public static double[][] InverseMatrix(double[][] matrix)
        {
            if (matrix.Length != 0 && matrix.Length != matrix[0].Length)
                throw new ArgumentException("Non-square matrix");

            if (Math.Abs(MatrixDeterminant(matrix)) < 0.000001)
                throw new ArgumentException("Degenerate matrix");

            int n = matrix[0].Length;
            double[][] reverseMatrix = new double[n][];

            for (int i = 0; i < n; ++i)
                reverseMatrix[i] = new double[n];

            for (int i = 0; i < n; ++i)
                reverseMatrix[i][i] = 1.0;

            for (int k = 0; k < n; ++k)
            {
                for (int i = 0; i < n; ++i)
                {
                    if (i != k)
                    {
                        for (int j = 0; j < n; ++j)
                        {
                            if (j != k)
                                matrix[i][j] -= matrix[k][j] * matrix[i][k] / matrix[k][k];
                            reverseMatrix[i][j] -= reverseMatrix[k][j] * matrix[i][k] / matrix[k][k];
                        }
                        matrix[i][k] = 0.0;
                    }
                }

                for (int i = 0; i < n; ++i)
                {
                    if (i != k)
                        matrix[k][i] /= matrix[k][k];
                    reverseMatrix[k][i] /= matrix[k][k];
                }
                matrix[k][k] = 1.0;
            }
            return reverseMatrix;
        }

        public static double[][] TransposeMatrix(double[][] matrix)
        {
            if (matrix.Length == 0)
                throw new ArgumentException("Empty matrix");

            double[][] result = new double[matrix[0].Length][];
            for (int k = 0; k < matrix[0].Length; ++k)
                result[k] = new double[matrix.Length];

            for (int i = 0; i < matrix.Length; ++i)
                for (int j = 0; j < matrix[0].Length; ++j)
                    result[j][i] = matrix[i][j];

            return result;
        }

        public static double[][] MultiplicateMatrix(double[][] a, double[][] b)
        {
            if (a[0].Length != b.Length)
                throw new Exception("Matrices can not be multiplied");

            double[][] product = new double[a.Length][];
            for (int i = 0; i < a.Length; ++i)
                product[i] = new double[b[0].Length];

            for (int i = 0; i < a.Length; ++i)
                for (int j = 0; j < b[0].Length; ++j)
                    for (int k = 0; k < b.Length; ++k)
                        product[i][j] += a[i][k] * b[k][j];

            return product;
        }

        public static double[] MultiplicateMatrix(double[] a, double[] b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors can not be multiplied");

            var product = new List<double>();
            for (int i = 0; i < a.Length; ++i)
                product.Add(a[i] * b[i]);
            return product.ToArray();
        }

        public static double[] MatrixVectorMultiplication(double[][] matrix, double[] vector)
        {
            if (matrix.Length != 0 && vector.Length != matrix[0].Length)
                throw new ArgumentException("Matrices can not be multiplied");

            double[] product = new double[matrix.Length];

            for (int i = 0; i < matrix.Length; ++i)
                for (int k = 0; k < matrix[i].Length; ++k)
                    product[i] += matrix[i][k] * vector[k];
            return product;
        }

        public static double SumOfVectorCoordinatesSquares(double[] vector)
        {
            double sum = 0;
            for (int i = 0; i < vector.Length; ++i)
                sum += vector[i] * vector[i];
            return sum;
        }

        public static double[] VectorsDifference(double[] a, double[] b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Lengths of the vectors are not equal");

            double[] difference = new double[a.Length];
            for (int i = 0; i < a.Length; ++i)
                difference[i] = a[i] - b[i];
            return difference;
        }

        public static double ScalarProductOfVectors(double[] a, double[] b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Lengths of the vectors are not equal");

            double product = 0;
            for (int i = 0; i < a.Length; i++)
                product += a[i] * b[i];
            return product;
        }

        public static double[] TransposeMatrixVectorProduct(double[] vector, double[][] matrix)
        {
            if (matrix.Length != 0 && (vector.Length != matrix[0].Length))
                throw new ArgumentException("Matrices can not be multiplied");

            double[] product = new double[vector.Length];
            for (int i = 0; i < matrix.Length; ++i)
                for (int k = 0; k < matrix[i].Length; ++k)
                    product[i] += vector[k] * matrix[i][k];
            return TransposeVector(product);
        }

        private static void DecomposeIntoLU(double[][] a, out double[][] l, out double[][] u)
        {
            l = new double[a.Length][];
            for (int i = 0; i < a.Length; ++i)
                l[i] = new double[a.Length];

            u = new double[a.Length][];
            for (int i = 0; i < a.Length; ++i)
                u[i] = new double[a.Length];

            int n = a.Length;
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    u[0][i] = a[0][i];
                    l[i][0] = a[i][0] / u[0][0];

                    double sum = .0;
                    for (int k = 0; k < i; ++k)
                    {
                        sum += l[i][k] * u[k][j];
                    }
                    u[i][j] = a[i][j] - sum;
                    if (i > j)
                    {
                        l[j][i] = 0;
                    }
                    else
                    {
                        sum = .0;
                        for (int k = 0; k < i; ++k)
                        {
                            sum += l[j][k] * u[k][i];
                        }
                        l[j][i] = (a[j][i] - sum) / u[i][i];
                    }
                }
            }
        }

        private static double UpperTriangularMatrixDeterminant(double[][] matrix)
        {
            double determinant = 1;
            for (int i = 0; i < matrix.Length; ++i)
                determinant *= matrix[i][i];
            return determinant;
        }

        private static double[][] DeleteRowAndColumnInMatrix(double[][] matrix, int row, int column)
        {
            int n = matrix.Length - 1;
            double[][] result = new double[n][];
            for (int i = 0; i < n; ++i)
                result[i] = new double[n];

            int x = 0;
            for (int i = 0; i < n + 1; ++i)
            {
                if (i != row)
                {
                    int y = 0;
                    for (int j = 0; j < n + 1; ++j)
                    {
                        if (j != column)
                        {
                            result[x][y] = matrix[i][j];
                            y++;
                        }
                    }
                    x++;
                }
            }
            return result;
        }
    }
}
