using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TVMS
{
    public class Forecast
    {
        public double[] Value { get; private set; }

        private double[][] dataMatrix;
        private double[] x0 = {1, 0.92, 0.83, 5.82, 1.2, 4.25, 1.01, 1.01, 1};

        public Forecast(double[][] matrix)
        {
            dataMatrix = matrix;
            Prognoz();
        }

        private void Prognoz()
        {
            int n = dataMatrix.Length, m = dataMatrix[0].Length;
            double[] result = new double[m];
            double[] Y = new double[n];
            double[][] A1 = new double[n][];
            for (int i = 0; i < n; i++)
            {
                A1[i] = new double[m];
                Y[i] = dataMatrix[i][0];
            }

            for (int i = 0; i < n; i++)
                A1[i][0] = 1;

            for (int i = 0; i < n; i++)
                for (int j = 1; j < m; j++)
                    A1[i][j] = dataMatrix[i][j];

            result = MatrixFunction.MatrixVectorMultiplication((MatrixFunction.MultiplicateMatrix(MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1)), MatrixFunction.TransposeMatrix(A1))), Y);
            double Qost = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.VectorsDifference(Y, MatrixFunction.MatrixVectorMultiplication(A1, result)));
            double s = Math.Sqrt(Qost / (n - m - 1));
            double t = 4.587;
            double[] a = new double[2];
            a[0] = MatrixFunction.ScalarProductOfVectors(x0, result) - t * s * Math.Sqrt(MatrixFunction.ScalarProductOfVectors(MatrixFunction.TransposeMatrixVectorProduct(x0, MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1))), x0));
            a[1] = MatrixFunction.ScalarProductOfVectors(x0, result) + t * s * Math.Sqrt(MatrixFunction.ScalarProductOfVectors(MatrixFunction.TransposeMatrixVectorProduct(x0, MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1))), x0));
            Value = a;
        }
    }
}
