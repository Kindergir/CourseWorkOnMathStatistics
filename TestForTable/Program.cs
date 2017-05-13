using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace TestForTable
{
    class Program
    {
        static void Main(string[] args)
        {
            int counter = 0;
            string line;

            StreamReader file = new StreamReader(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata.csv");
            StreamWriter sw = new StreamWriter(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata_1.csv");


            double a = 0, b, c, d, e, f, g, h, t, y, u;
            double mina = 4230326280, minb = 4230326280, minc = 4230326280, mind = 4230326280, mine = 4230326280, minf = 4230326280, ming = 4230326280,
                minh = 4230326280, mint = 4230326280, miny = 4230326280, minu = 4230326280;
            double maxa = -1, maxb = -1, maxc = -1, maxd = -1, maxe = -1, maxf = -1, maxg = -1, maxh = -1, maxt = -1, maxy = -1, maxu = -1;

            while ((line = file.ReadLine()) != null)
            {
                var cur = line.Split(',');

                //if (counter == 0)
                //{
                //    sw.WriteLine("Title;Genres;Duration;People participated in the shooting;Director score;Actor1 score;Actor2 score;Actor3 score;People waiting;Requests on network;Budget;Gross;IMDB-score");
                //}
                //else if (cur[11] != "" && cur[9] != "" && cur[3] != "" && cur[2] != "" && cur[4] != "" && cur[7] != "" &&
                //     cur[24] != "" && cur[5] != "" && cur[13] != "" && cur[12] != "" && cur[22] != "" && cur[8] != "" && cur[25] != "")
                //sw.WriteLine(cur[11] + ';' + cur[9] + ';' + cur[3] + ';' + cur[2] + ';' + cur[4] + ';' + cur[7] + ';'
                //    + cur[24] + ';' + cur[5] + ';' + cur[13] + ';' + cur[12] + ';' + cur[22] + ';' + cur[8] + ';' + cur[25].Replace('.', ','));

                if (counter == 0)
                {
                    //sw.WriteLine("Duration;People participated in the shooting;Director score;Actor1 score;Actor2 score;Actor3 score;People waiting;Requests on network;Budget;Gross;IMDB-score");
                }
                else if (cur[3] != "" && cur[2] != "" && cur[4] != "" && cur[7] != "" &&
                     cur[24] != "" && cur[5] != "" && cur[13] != "" && cur[12] != "" && cur[22] != "" && cur[8] != "" && cur[25] != ""
                    && Double.TryParse(cur[3], out a) && Double.TryParse(cur[2], out b) && Double.TryParse(cur[4], out c) && Double.TryParse(cur[7], out d) &&
                     Double.TryParse(cur[24], out e) && Double.TryParse(cur[5], out f) && Double.TryParse(cur[13], out g) && Double.TryParse(cur[12], out h)
                    && Double.TryParse(cur[22], out y) && Double.TryParse(cur[8], out u) && Double.TryParse(cur[25].Replace('.', ','), out t))
                {
                    if (a < mina)
                        mina = a;
                    if (b < minb)
                        minb = b;
                    if (c < minc)
                        minc = c;
                    if (d < mind)
                        mind = d;
                    if (e < mine)
                        mine = e;
                    if (f < minf)
                        minf = f;
                    if (g < ming)
                        ming = g;
                    if (h < minh)
                        minh = h;
                    if (y < miny)
                        miny = y;
                    if (u < minu)
                        minu = u;
                    if (t < mint)
                        mint = t;
                    if (a > maxa)
                        maxa = a;
                    if (b > maxb)
                        maxb = b;
                    if (c > maxc)
                        maxc = c;
                    if (d > maxd)
                        maxd = d;
                    if (e > maxe)
                        maxe = e;
                    if (f > maxf)
                        maxf = f;
                    if (g > maxg)
                        maxg = g;
                    if (h > maxh)
                        maxh = h;
                    if (y > maxy)
                        maxy = y;
                    if (u > maxu)
                        maxu = u;
                    if (t > maxt)
                        maxt = t;
                }
                counter++;
            }
            file.Close();

            //StreamReader file = new StreamReader(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata.csv");
            //StreamWriter sw = new StreamWriter(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata_1.csv");
            StreamReader file1 = new StreamReader(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata.csv");
            while ((line = file1.ReadLine()) != null)
            {
                var cur = line.Split(',');

                //if (counter == 0)
                //{
                //    sw.WriteLine("Title;Genres;Duration;People participated in the shooting;Director score;Actor1 score;Actor2 score;Actor3 score;People waiting;Requests on network;Budget;Gross;IMDB-score");
                //}
                //else if (cur[11] != "" && cur[9] != "" && cur[3] != "" && cur[2] != "" && cur[4] != "" && cur[7] != "" &&
                //     cur[24] != "" && cur[5] != "" && cur[13] != "" && cur[12] != "" && cur[22] != "" && cur[8] != "" && cur[25] != "")
                //sw.WriteLine(cur[11] + ';' + cur[9] + ';' + cur[3] + ';' + cur[2] + ';' + cur[4] + ';' + cur[7] + ';'
                //    + cur[24] + ';' + cur[5] + ';' + cur[13] + ';' + cur[12] + ';' + cur[22] + ';' + cur[8] + ';' + cur[25].Replace('.', ','));
                //double a = 0, b, c, d, e, f, g, h, t, y, u;
                //double mina = 4230326280, minb = 4230326280, minc = 4230326280, mind = 4230326280, mine = 4230326280, minf = 4230326280, ming = 4230326280,
                //    minh = 4230326280, mint = 4230326280, miny = 4230326280, minu = 4230326280;
                //double maxa = -1, maxb = -1, maxc = -1, maxd = -1, maxe = -1, maxf = -1, maxg = -1, maxh = -1, maxt = -1, maxy = -1, maxu = -1;
                if (counter == 0)
                {
                    //sw.WriteLine("Duration;People participated in the shooting;Director score;Actor1 score;Actor2 score;Actor3 score;People waiting;Requests on network;Budget;Gross;IMDB-score");
                }
                else if (cur[3] != "" && cur[2] != "" && cur[4] != "" && cur[7] != "" &&
                     cur[24] != "" && cur[5] != "" && cur[13] != "" && cur[12] != "" && cur[22] != "" && cur[8] != "" && cur[25] != ""
                    && Double.TryParse(cur[3], out a) && Double.TryParse(cur[2], out b) && Double.TryParse(cur[4], out c) && Double.TryParse(cur[7], out d) &&
                     Double.TryParse(cur[24], out e) && Double.TryParse(cur[5], out f) && Double.TryParse(cur[13], out g) && Double.TryParse(cur[12], out h)
                    && Double.TryParse(cur[22], out y) && Double.TryParse(cur[8], out u) && Double.TryParse(cur[25].Replace('.', ','), out t)
                    && cur[3] != "0" && cur[2] != "0" && cur[4] != "0" && cur[7] != "0" &&
                     cur[24] != "0" && cur[5] != "0" && cur[13] != "0" && cur[12] != "0" && cur[22] != "0" && cur[8] != "0" && cur[25] != "0")
                //sw.WriteLine(cur[3] + ';' + cur[2] + ';' + cur[4] + ';' + cur[7] + ';'
                //    + cur[24] + ';' + cur[5] + ';' + cur[13] + ';' + cur[12] + ';' + cur[22] + ';' + cur[8] + ';' + cur[25].Replace('.', ','));
                    sw.WriteLine((a / (maxa - mina)).ToString() + ';' + (b / (maxb - minb)).ToString() + ';' + (c / (maxc - minc)).ToString() + ';' + 
                        (d / (maxd - mind)).ToString() + ';'
                    + (e / (maxe - mine)).ToString() + ';' + (f / (maxf - minf)).ToString() + ';' + (g / (maxg - ming)).ToString() + ';' + (h / (maxh - minh)).ToString() 
                    + ';' + (y / (maxy - miny)).ToString() + ';'
                    + (u / (maxu - minu)).ToString() + ';' + (t / (maxt - mint)).ToString());
                counter++;
            }
            file.Close();
        } 
    }
}
