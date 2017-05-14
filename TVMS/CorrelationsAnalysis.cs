using System;
using System.Linq;

namespace TVMS
{
    public class CorrelationsAnalysis
    {
        public double[][] PairCorrelationsMatrix { get; private set; }
        public double[][] MatrixSignificanceFactors { get; private set; }
        public double[][] PartialCorrelationsMatrix { get; private set; }
        public double SelectiveMultipleCoefficient { get; private set; }
        public double DeterminationCoefficient { get; private set; }
        public double[] ParametersSignificanceFactors { get; private set; }

        private double[][] dataMatrix;

        public CorrelationsAnalysis(double[][] matrix, int resultParameterNumber)
        {
            dataMatrix = matrix;
            CalcPairCorrelationsMatrix();
            CalcSignificanceFactors();
            CalcPartialCorrelationsMatrix();
            CalcSelectiveMultipleCoefficient(resultParameterNumber);
            CalcDeterminationCoefficient();
            CalcParametersSignificanceFactors();
        }

        private void CalcPairCorrelationsMatrix()
        {
            var averageValues = AverageValuesForAllCoefficients(dataMatrix);
            var squaredDevivations = SumSquaredDevivatiosForAllCofficients(dataMatrix);
            var devivations = SumDevivatiosForAllCofficients(dataMatrix);

            int coefCount = dataMatrix.Length;
            PairCorrelationsMatrix = new double[coefCount][];
            for (int k = 0; k < coefCount; k++)
                PairCorrelationsMatrix[k] = new double[coefCount];

            for (int i = 0; i < coefCount; ++i)
            {
                for (int j = i; j < coefCount; ++j)
                {
                    if (i == j)
                        PairCorrelationsMatrix[i][j] = 1;
                    else
                    {
                        var xAverage = averageValues[i];
                        var yAverage = averageValues[j];
                        var xDevivationsSum = devivations[i];
                        var yDevivationsSum = devivations[j];
                        var xSquaredDevivationsSum = squaredDevivations[i];
                        var ySquaredDevivationsSum = squaredDevivations[j];

                        PairCorrelationsMatrix[i][j] = (xDevivationsSum - xAverage) * (yDevivationsSum - yAverage)
                                / Math.Sqrt(xSquaredDevivationsSum * ySquaredDevivationsSum);
                    }
                }
            }
        }

        private void CalcSignificanceFactors()
        {
            MatrixSignificanceFactors = new double[PairCorrelationsMatrix.Length][];
            for (int i = 0; i < PairCorrelationsMatrix.Length; ++i)
                MatrixSignificanceFactors[i] = new double[PairCorrelationsMatrix.Length];

            for (int i = 0; i < PairCorrelationsMatrix.Length; ++i)
            {
                for (int j = 0; j < PairCorrelationsMatrix.Length; ++j)
                {
                    if (i == j)
                        MatrixSignificanceFactors[i][j] = 1;
                    else
                        MatrixSignificanceFactors[i][j] = Math.Round(Math.Sqrt(Math.Pow(PairCorrelationsMatrix[i][j], 2)*(dataMatrix[i].Length - 2))/
                                (1 - Math.Pow(PairCorrelationsMatrix[i][j], 2)), 3);
                }
            }
        }

        private void CalcPartialCorrelationsMatrix()
        {
            int n = PairCorrelationsMatrix.Length;
            PartialCorrelationsMatrix = new double[n][];
            for (int i = 0; i < n; i++)
                PartialCorrelationsMatrix[i] = new double[n];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    PartialCorrelationsMatrix[i][j] = MatrixFunction.AlgebraicComplement(PairCorrelationsMatrix, i, j)
                        / Math.Sqrt(MatrixFunction.AlgebraicComplement(PairCorrelationsMatrix, i, i) 
                        * MatrixFunction.AlgebraicComplement(PairCorrelationsMatrix, j, j));
        }

        private void CalcSelectiveMultipleCoefficient(int resultParameterNumber)
        {
            SelectiveMultipleCoefficient = Math.Sqrt(1 - MatrixFunction.MatrixDeterminant(PairCorrelationsMatrix)
                / MatrixFunction.AlgebraicComplement(PairCorrelationsMatrix, resultParameterNumber, resultParameterNumber));
        }

        private void CalcDeterminationCoefficient()
        {
            DeterminationCoefficient = Math.Pow(SelectiveMultipleCoefficient, 2);
        }

        private void CalcParametersSignificanceFactors()
        {
            int parametersCount = PairCorrelationsMatrix.Length;
            ParametersSignificanceFactors = new double[parametersCount];
            for (int i = 0; i < parametersCount; ++i)
                ParametersSignificanceFactors[i] = Math.Pow(SelectiveMultipleCoefficient, 2)
                    / ((1 - Math.Pow(SelectiveMultipleCoefficient, 2)) / (parametersCount - 2));
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
