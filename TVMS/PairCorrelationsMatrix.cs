using System;
using System.Linq;

namespace TVMS
{
    public class PairCorrelationsMatrix
    {
        public double[][] CorrelationsMatrix { get; private set; }
        public double[][] SignificanceFactorsMatrix { get; private set; }

        private double[][] dataMatrix;

        public PairCorrelationsMatrix(double[][] matrix)
        {
            dataMatrix = matrix;
            CalcCorrelationsMatrix();
            CalcSignificanceFactors();
        }

        public void CalcCorrelationsMatrix()
        {
            var averageValues = AverageValuesForAllCoefficients(dataMatrix);
            var squaredDevivations = SumSquaredDevivatiosForAllCofficients(dataMatrix);
            var devivations = SumDevivatiosForAllCofficients(dataMatrix);

            int coefCount = dataMatrix.Length;
            CorrelationsMatrix = new double[coefCount][];
            for (int k = 0; k < coefCount; k++)
                CorrelationsMatrix[k] = new double[coefCount];

            for (int i = 0; i < coefCount; ++i)
            {
                for (int j = i; j < coefCount; ++j)
                {
                    if (i == j)
                        CorrelationsMatrix[i][j] = 1;
                    else
                    {
                        var xAverage = averageValues[i];
                        var yAverage = averageValues[j];
                        var xDevivationsSum = devivations[i];
                        var yDevivationsSum = devivations[j];
                        var xSquaredDevivationsSum = squaredDevivations[i];
                        var ySquaredDevivationsSum = squaredDevivations[j];

                        CorrelationsMatrix[i][j] = (xDevivationsSum - xAverage) * (yDevivationsSum - yAverage)
                                / Math.Sqrt(xSquaredDevivationsSum * ySquaredDevivationsSum);
                    }
                }
            }
        }

        private void CalcSignificanceFactors()
        {
            SignificanceFactorsMatrix = new double[CorrelationsMatrix.Length][];
            for (int i = 0; i < CorrelationsMatrix.Length; ++i)
                SignificanceFactorsMatrix[i] = new double[CorrelationsMatrix.Length];

            for (int i = 0; i < CorrelationsMatrix.Length; ++i)
            {
                for (int j = 0; j < CorrelationsMatrix.Length; ++j)
                {
                    if (i == j)
                        SignificanceFactorsMatrix[i][j] = 1;
                    else
                        SignificanceFactorsMatrix[i][j] = Math.Round(Math.Sqrt(Math.Pow(CorrelationsMatrix[i][j], 2)*(dataMatrix[i].Length - 2))/
                                (1 - Math.Pow(CorrelationsMatrix[i][j], 2)), 3);
                }
            }
        }

        private double SumSquaredDeviationsOfQuantity(double[] values)
        {
            double x = 0;
            var averageValue = values.Average();
            for (int i = 0; i < values.Length; ++i)
                x += Math.Pow(values[i] - averageValue, 2);
            return x;
        }

        private double SumDeviationsOfQuantity(double[] values)
        {
            double x = 0;
            var averageValue = values.Average();
            for (int i = 0; i < values.Length; ++i)
                x += values[i] - averageValue;
            return x;
        }

        private double[] AverageValuesForAllCoefficients(double[][] matrix)
        {
            double[] values = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                values[i] = matrix[i].Average();
            return values;
        }

        private double[] SumDevivatiosForAllCofficients(double[][] matrix)
        {
            double[] devivations = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                devivations[i] = SumDeviationsOfQuantity(matrix[i]);
            return devivations;
        }

        private double[] SumSquaredDevivatiosForAllCofficients(double[][] matrix)
        {
            double[] devivations = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                devivations[i] = SumSquaredDeviationsOfQuantity(matrix[i]);
            return devivations;
        }
    }
}
