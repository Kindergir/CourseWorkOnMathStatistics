using System;
using System.Collections.Generic;
using System.Linq;

namespace TVMS
{
    public class PearsonConsentCriterion
    {
        public double[] PearsonCriterionValue { get; }
        public double[] AverageValueX { get; }
        public double[][] Intervals { get; }
        public int[][] HitsInIntervalsCount { get; }
        public double[][] HitsInIntervalsProbability { get; }
        public double[] MeanSquareDeviation { get; }

        public PearsonConsentCriterion(double[][] dataMatrix, Dictionary<double, double> laplasMatrix)
        {
            this.laplasMatrix = laplasMatrix;
            int parametersCount = dataMatrix.Length;
            Intervals = new double[parametersCount][];
            HitsInIntervalsCount = new int[parametersCount][];
            MeanSquareDeviation = new double[parametersCount];
            HitsInIntervalsProbability = new double[parametersCount][];
            AverageValueX = new double[parametersCount];
            PearsonCriterionValue = new double[parametersCount];

            for (int i = 0; i < parametersCount; ++i)
            {
                Intervals[i] = BuildItervals(dataMatrix[i], parametersCount).ToArray();
                HitsInIntervalsCount[i] = CalcHitsInIntervalsCount(dataMatrix[i], Intervals[i]);
                AverageValueX[i] = CalcAverageX(Intervals[i].ToList(), HitsInIntervalsCount[i]);
                MeanSquareDeviation[i] = CalcMeanSquareDeviation(Intervals[i], HitsInIntervalsCount[i], AverageValueX[i], dataMatrix[0].Length);
                HitsInIntervalsProbability[i] = CalcProbabilityHitsInInterval(Intervals[i], AverageValueX[i], MeanSquareDeviation[i], laplasMatrix);
                PearsonCriterionValue[i] = PearsonTestForOneParameter(HitsInIntervalsCount[i], HitsInIntervalsProbability[i], Intervals[i].Length);
            }
        }

        private Dictionary<double, double> laplasMatrix;

        private double PearsonTestForOneParameter(int[] hitsCount, double[] hitsProbability, int intervalsCount)
        {
            double sum = 0;
            for (int i = 0; i < hitsCount.Length; ++i)
                sum += Math.Pow(hitsCount[i] - intervalsCount * hitsProbability[i], 2.0)
                    / (intervalsCount * hitsProbability[i]);
            return sum;
        }

        private List<double> BuildItervals(double[] d, int varCount)
        {
            var length = SpotLengthInterval(d, varCount);
            List<double> a = new List<double>();
            double a0 = d.Min() - length / 2;
            double bound = a0 + length;
            a.Add(d.Min() - length / 2);
            a.Add(bound);
            while (bound <= d.Max())
            {
                bound += length;
                a.Add(bound);
            }
            return a;
        }

        private double SpotLengthInterval(double[] d, int n)
        {
            return (d.Max() - d.Min()) / (1 + 3.321 * Math.Log10(n));
        }

        private double CalcAverageX(List<double> a, int[] n)
        {
            double[] x = new double[a.Count - 1];
            for (int j = 0; j < a.Count - 1; ++j)
                x[j] = (a[j] + a[j + 1]) / 2;

            double X = 0, F = 0;
            for (int i = 0; i < n.Length; ++i)
            {
                X += x[i] * n[i];
                F += n[i];
            }
            return X / F;
        }

        private double CalcMeanSquareDeviation(double[] a, int[] n, double x, int m)
        {
            double s = 0;
            double[] X = new double[a.Length - 1];
            for (int j = 0; j < a.Length - 1; j++)
            {
                X[j] = (a[j] + a[j + 1]) / 2;
            }

            for (int i = 0; i < n.Length; i++)
            {
                s += Math.Pow(X[i] - x, 2) * n[i];
            }

            return Math.Sqrt(s / m);
        }

        private int[] CalcHitsInIntervalsCount(double[] values, double[] intervals)
        {
            int[] n = new int[intervals.Length - 1];

            for (int i = 0; i < intervals.Length - 1; ++i)
                for (int j = 0; j < values.Length; ++j)
                    if (values[j] > intervals[i] && values[j] <= intervals[i + 1])
                        n[i]++;
            return n;
        }

        private double[] CalcProbabilityHitsInInterval(double[] a, double X, double s, Dictionary<double, double> L)
        {
            double[] Ni = new double[a.Length - 1];
            double t1, t2;
            for (int i = 0; i < a.Length - 1; i++)
            {
                t1 = (a[i + 1] - X) / s;
                t2 = (a[i] - X) / s;
                if (t1 < 0)
                {
                    L.TryGetValue(Math.Round(-t1, 2), out t1);
                    t1 *= -1;
                }
                else
                {
                    L.TryGetValue(Math.Round(t1, 2), out t1);
                }

                if (t2 < 0)
                {
                    L.TryGetValue(Math.Round(-t2, 2), out t2);
                    t2 *= -1;
                }
                else
                {
                    L.TryGetValue(Math.Round(t2, 2), out t2);
                }
                Ni[i] = t1 - t2;
            }
            return Ni;
        }
    }
}
