using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TVMS
{
    public static class PartialCorrelationsMatrix
    {
        public static double[][] GetForMatrix(double[][] pairCorrelations)
        {
            int n = pairCorrelations.Length;
            double[][] partialCorrelations = new double[n][];
            for (int i = 0; i < n; i++)
                partialCorrelations[i] = new double[n];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    partialCorrelations[i][j] = MatrixFunction.AlgebraicComplement(pairCorrelations, i, j)
                        / Math.Sqrt(MatrixFunction.AlgebraicComplement(pairCorrelations, i, i) * MatrixFunction.AlgebraicComplement(pairCorrelations, j, j));
            return partialCorrelations;
        }
    }
}
