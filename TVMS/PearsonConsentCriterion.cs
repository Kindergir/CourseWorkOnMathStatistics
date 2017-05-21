using System;
using System.Collections.Generic;
using System.Linq;

namespace TVMS
{
    public class PearsonConsentCriterion
    {
        public double PearsonCriterionValue { get; private set; }
        public double AverageValueX { get; private set; }
        public double[] Intervals { get; private set; }
        public int[] HitsInIntervalsCount { get; private set; }
        public double[] HitsInIntervalsProbability { get; private set; }
        public double MeanSquareDeviation { get; private set; }

        private Dictionary<double, double> laplasMatrix;

        public PearsonConsentCriterion(double[] parameterValues, Dictionary<double, double> laplasMatrix)
        {
            this.laplasMatrix = laplasMatrix;
            BuildItervals(parameterValues);
            CalcHitsInIntervalsCount(parameterValues);
            CalcAverageX();
            CalcMeanSquareDeviation(parameterValues.Length);
            CalcProbabilityHitsInInterval();
            PearsonTestForOneParameter(parameterValues.Length);
        }

        private void PearsonTestForOneParameter(int measureCount)
        {
            double sum = 0;
            for (int i = 0; i < HitsInIntervalsCount.Length; ++i)
                if (HitsInIntervalsProbability[i] != 0)
                    sum += Math.Pow(HitsInIntervalsCount[i] - measureCount * HitsInIntervalsProbability[i], 2.0)
                        / (measureCount * HitsInIntervalsProbability[i]);
            PearsonCriterionValue = sum;
        }

        private void BuildItervals(double[] parameterValues)
        {
            var length = SpotLengthInterval(parameterValues);
            List<double> a = new List<double>();
            double a0 = parameterValues.Min() - length / 2;
            double bound = a0 + length;
            a.Add(parameterValues.Min() - length / 2);
            a.Add(bound);
            while (bound <= parameterValues.Max())
            {
                bound += length;
                a.Add(bound);
            }
            Intervals = a.ToArray();
        }

        private double SpotLengthInterval(double[] parameterValues)
        {
            return (parameterValues.Max() - parameterValues.Min()) / (1 + 3.321 * Math.Log10(parameterValues.Length));
        }

        private void CalcAverageX()
        {
            double[] midpoints = new double[Intervals.Length - 1];
            for (int j = 0; j < Intervals.Length - 1; ++j)
                midpoints[j] = (Intervals[j] + Intervals[j + 1]) / 2;

            double meanValuesSum = 0, F = 0;
            for (int i = 0; i < HitsInIntervalsCount.Length; ++i)
            {
                meanValuesSum += midpoints[i] * HitsInIntervalsCount[i];
                F += HitsInIntervalsCount[i];
            }
            AverageValueX = meanValuesSum / F;
        }

        private void CalcMeanSquareDeviation(int measureCount)
        {
            double s = 0;
            double[] midpoints = new double[Intervals.Length - 1];
            for (int j = 0; j < Intervals.Length - 1; j++)
                midpoints[j] = (Intervals[j] + Intervals[j + 1]) / 2;

            for (int i = 0; i < HitsInIntervalsCount.Length; i++)
                s += Math.Pow(midpoints[i] - AverageValueX, 2) * HitsInIntervalsCount[i];
            MeanSquareDeviation = Math.Sqrt(s / measureCount);
        }

        private void CalcHitsInIntervalsCount(double[] values)
        {
            List<int> hits = new List<int>();
            List<double> intervals = new List<double>();
            intervals.Add(Intervals[0]);

            int add = 0;
            for (int i = 0; i < Intervals.Length - 1; ++i)
            {
                var cnt = values.Count(x => x > Intervals[i] && x <= Intervals[i + 1]);
                if (cnt != 0)
                {
                    if (cnt >= 5)
                    {
                        hits.Add(cnt + add);
                        intervals.Add(Intervals[i + 1]);
                        add = 0;
                    }
                    else
                        add = cnt;
                }
            }
            if (add == 0)
            {
                HitsInIntervalsCount = hits.ToArray();
                if (intervals[intervals.Count - 1] != Intervals[Intervals.Length - 1])
                    intervals.Add(Intervals[Intervals.Length - 1]);
            }
            else
            {
                HitsInIntervalsCount = hits.ToArray();
                intervals[intervals.Count - 1] = Intervals[Intervals.Length - 1];
            }
            Intervals = intervals.ToArray();
        }

        private void CalcProbabilityHitsInInterval()
        {
            HitsInIntervalsProbability = new double[Intervals.Length - 1];
            double t1, t2;
            for (int i = 0; i < Intervals.Length - 1; i++)
            {
                t1 = Intervals[i + 1] - AverageValueX/MeanSquareDeviation;
                t2 = Intervals[i] - AverageValueX/MeanSquareDeviation;

                if (t1 < 0)
                {
                    if (laplasMatrix.ContainsKey(Math.Round(-1 * t1, 2)))
                        t1 *= laplasMatrix[Math.Round(-1 * t1, 2)] * -1;
                    else
                        t1 *= -0.5;
                }
                else
                {
                    if (laplasMatrix.ContainsKey(Math.Round(t1, 2)))
                        t1 *= laplasMatrix[Math.Round(t1, 2)];
                    else
                        t1 *= 0.5;
                }

                if (t2 < 0)
                {
                    if (laplasMatrix.ContainsKey(Math.Round(-1 * t2, 2)))
                        t2 *= laplasMatrix[Math.Round(-1 * t2, 2)] * -1;
                    else
                        t2 *= -0.5;
                }
                else
                {
                    if (laplasMatrix.ContainsKey(Math.Round(t2, 2)))
                        t2 *= laplasMatrix[Math.Round(t2, 2)];
                    else
                        t2 *= 0.5;
                }

                HitsInIntervalsProbability[i] = Math.Abs(t1 - t2);
            }
        }
    }
}
