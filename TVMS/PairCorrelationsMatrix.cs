using System;
using System.Linq;

namespace TVMS
{
    public class PairCorrelationsMatrix
    {
        public static double[][] GetForMatrix(double[][] matrix)
        {
            var averageValues = AverageValuesForAllCoefficients(matrix);
            var squaredDevivations = SumSquaredDevivatiosForAllCofficients(matrix);
            var devivations = SumDevivatiosForAllCofficients(matrix);

            int coefCount = matrix.Length;
            double[][] pairCorrelationsMatrix = new double[coefCount][];
            for (int k = 0; k < coefCount; k++)
                pairCorrelationsMatrix[k] = new double[coefCount];

            for (int i = 0; i < coefCount; ++i)
            {
                for (int j = i; j < coefCount; ++j)
                {
                    if (i == j)
                        pairCorrelationsMatrix[i][j] = 1;
                    else
                    {
                        var xAverage = averageValues[i];
                        var yAverage = averageValues[j];
                        var xDevivationsSum = devivations[i];
                        var yDevivationsSum = devivations[j];
                        var xSquaredDevivationsSum = squaredDevivations[i];
                        var ySquaredDevivationsSum = squaredDevivations[j];

                        pairCorrelationsMatrix[i][j] = (xDevivationsSum - xAverage) * (yDevivationsSum - yAverage)
                                / Math.Sqrt(xSquaredDevivationsSum * ySquaredDevivationsSum);
                    }
                }
            }
            return pairCorrelationsMatrix;
        }

        private static double SumSquaredDeviationsOfQuantity(double[] values)
        {
            double x = 0;
            var averageValue = values.Average();
            for (int i = 0; i < values.Length; ++i)
                x += Math.Pow(values[i] - averageValue, 2);
            return x;
        }

        private static double SumDeviationsOfQuantity(double[] values)
        {
            double x = 0;
            var averageValue = values.Average();
            for (int i = 0; i < values.Length; ++i)
                x += values[i] - averageValue;
            return x;
        }

        private static double[] AverageValuesForAllCoefficients(double[][] matrix)
        {
            double[] values = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                values[i] = matrix[i].Average();
            return values;
        }

        private static double[] SumDevivatiosForAllCofficients(double[][] matrix)
        {
            double[] devivations = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                devivations[i] = SumDeviationsOfQuantity(matrix[i]);
            return devivations;
        }

        private static double[] SumSquaredDevivatiosForAllCofficients(double[][] matrix)
        {
            double[] devivations = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                devivations[i] = SumSquaredDeviationsOfQuantity(matrix[i]);
            return devivations;
        }
    }
}
