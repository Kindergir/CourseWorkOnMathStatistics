using System;
using System.Collections.Generic;
using System.Linq;

namespace MathStatisticsLibrary
{
    public class Forecast
    {
        public double[] Value { get; private set; }

        private double[][] dataMatrix;
        private double[] x0 = {1, 0.92, 0.83, 5.82, 1.2, 4.25, 1.01, 1.01, 1, 1};
        private int resultParameterNumber;

        public Forecast(double[][] matrix, int resultParameterNumber)
        {
            dataMatrix = matrix;
            this.resultParameterNumber = resultParameterNumber;
            Prognoz();
        }

        private void Prognoz()
        {
            var transpMatrix = MatrixFunction.TransposeMatrix(dataMatrix);
            int parametersCount = dataMatrix[0].Length, measuresCount = dataMatrix.Length;

            double[] resultParameter = transpMatrix[resultParameterNumber];
            List<double>[] withoutResultParameter = new List<double>[measuresCount];

            for (int i = 0; i < measuresCount; i++)
            {
                withoutResultParameter[i] = new List<double>();
                withoutResultParameter[i].Add(1);

                for (int j = 0; j < parametersCount; ++j)
                    if (j != resultParameterNumber)
                        withoutResultParameter[i].Add(dataMatrix[i][j]);
            }

            var xMatrix = withoutResultParameter.Select(o => o.ToArray()).ToArray();
            var transpWithoutMatrix = MatrixFunction.TransposeMatrix(xMatrix);
            RegressionAnalysis ra = new RegressionAnalysis(dataMatrix, resultParameterNumber);

            var result = ra.RegressionCoefficients;
            double s = ra.ResidualDispersion;
            double t = 4.587;
            double[] a = new double[2];
            a[0] = MatrixFunction.ScalarProductOfVectors(x0, result) - t * s * Math.Sqrt(MatrixFunction.ScalarProductOfVectors(MatrixFunction.TransposeMatrixVectorProduct(x0, MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(transpWithoutMatrix, xMatrix))), x0));
            a[1] = MatrixFunction.ScalarProductOfVectors(x0, result) + t * s * Math.Sqrt(MatrixFunction.ScalarProductOfVectors(MatrixFunction.TransposeMatrixVectorProduct(x0, MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(transpWithoutMatrix, xMatrix))), x0));
            Value = a;
        }

    }
}
