using System;
using System.Collections.Generic;
using System.Linq;

namespace TVMS
{
    public static class ConsentCriterion
    {
        static List<double> PearsonChiSquaredTest(double[][] matrix, Dictionary<double, double> laplasMatrix)
        {
            var criteries = new List<double>();
            for (int i = 0; i < matrix.Length; ++i)
            {
                List<double> intervals = BuildItervals(matrix[i], matrix.Length);
                int intervalsCount = intervals.Count;
                int[] hitsInIntervalsCount = ValueHitsInIntervalsCount(matrix[i], intervals);
                double meanSquareDeviation = MeanSquareDeviation(intervals, hitsInIntervalsCount, Xsr(intervals, hitsInIntervalsCount), n);
                double[] hitsProbability = ni(intervals, Xsr(intervals, hitsInIntervalsCount), meanSquareDeviation, laplasMatrix);

                double sum = 0;
                for (int j = 0; j < hitsInIntervalsCount.Length; ++j)
                    sum += Math.Pow(hitsInIntervalsCount[i] - intervalsCount * hitsProbability[i], 2.0)
                           / (intervalsCount * hitsProbability[i]);
                criteries.Add(sum);
            }
            return criteries;
        }

        static List<double> BuildItervals(double[] d, int varCount)
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

        static double SpotLengthInterval(double[] d, int n)
        {
            return (d.Max() - d.Min()) / (1 + 3.321 * Math.Log10(n));
        }

        static double Xsr(List<double> a, int[] n)
        {
            double[] x = new double[a.Count - 1];
            for (int j = 0; j < a.Count - 1; j++)
            {
                x[j] = (a[j] + a[j + 1]) / 2;
            }

            double X = 0, F = 0;
            for (int i = 0; i < n.Length; i++)
            {
                X += x[i] * n[i];
                F += n[i];
            }
            return X / F;
        }

        static double MeanSquareDeviation(List<double> a, int[] n, double x, int m)
        {
            double s = 0;
            double[] X = new double[a.Count - 1];
            for (int j = 0; j < a.Count - 1; j++)
            {
                X[j] = (a[j] + a[j + 1]) / 2;
            }

            for (int i = 0; i < n.Length; i++)
            {
                s += Math.Pow(X[i] - x, 2) * n[i];
            }

            return Math.Sqrt(s / m);
        }

        static int[] ValueHitsInIntervalsCount(double[] values, List<double> intervals)
        {
            int[] n = new int[intervals.Count - 1];

            for (int i = 0; i < intervals.Count - 1; ++i)
                for (int j = 0; j < values.Length; ++j)
                    if (values[j] > intervals[i] && values[j] <= intervals[i + 1])
                        n[i]++;
            return n;
        }

        static double[] ni(List<double> a, double X, double s, Dictionary<double, double> L)
        {
            double[] Ni = new double[a.Count - 1];
            double t1, t2;
            for (int i = 0; i < a.Count - 1; i++)
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
