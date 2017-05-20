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
            var squaredDevivations = CalcSquaredDevivatiosForAllCofficients(dataMatrix);
            var devivations = CalcSumDevivatiosForAllCofficients(dataMatrix);

            int coefCount = dataMatrix.Length;
            int measurementsCount = dataMatrix[0].Length;
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
                        var xDevivationsSum = devivations[i];
                        var yDevivationsSum = devivations[j];
                        var xSquaredDevivation = squaredDevivations[i];
                        var ySquaredDevivation = squaredDevivations[j];

                        PairCorrelationsMatrix[i][j] = xDevivationsSum * yDevivationsSum
                                / xSquaredDevivation * ySquaredDevivation * measurementsCount;
                        PairCorrelationsMatrix[j][i] = PairCorrelationsMatrix[i][j];
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
                ParametersSignificanceFactors[i] = DeterminationCoefficient * Math.Sqrt(parametersCount - 2)
                    / 1 - DeterminationCoefficient;
        }

        private double[] CalcSumDevivatiosForAllCofficients(double[][] matrix)
        {
            double[] devivations = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                devivations[i] = SumDeviationsOfQuantity(matrix[i]);
            return devivations;
        }

        private double SumDeviationsOfQuantity(double[] values)
        {
            double x = 0;
            var averageValue = values.Average();
            for (int i = 0; i < values.Length; ++i)
                x += values[i] - averageValue;
            return x;
        }

        private double[] CalcSquaredDevivatiosForAllCofficients(double[][] matrix)
        {
            double[] devivations = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; ++i)
                devivations[i] = CalcQuantitySquaredDeviation(matrix[i]);
            return devivations;
        }

        private double CalcQuantitySquaredDeviation(double[] values)
        {
            double x = 0;
            var averageValue = values.Average();
            for (int i = 0; i < values.Length; ++i)
                x += Math.Pow(values[i] - averageValue, 2);
            return Math.Sqrt(1.0 / values.Length * x);
        }
    }
}
