#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//Un pixel este un triplet RGB
typedef struct
{
    unsigned char R, G, B;
} pixel;

void dimensiuni_imagine (char *nume_fisier, unsigned int *H, unsigned int *W, unsigned int *padding)
{
    //lungime si inaltime
    FILE *fin;
    fin = fopen (nume_fisier, "rb");
    if (fin == NULL)
    {
        printf ("Citire gresita.");
        return;
    }
    *H = 0;
    *W = 0;
    int i;
    unsigned int x = 0;
    for (i = 0; i <= 54; i = i + 2)
    {
        fread (&x, sizeof (int) / 2, 1, fin);
        if (i == 18)
            *W = x;
        if (i == 22)
            *H = x;
    }
    if ( (*W) % 4 != 0)
        *padding = 4 - (3 * (*W) ) % 4;
    else
        *padding = 0;
    fclose (fin);
}

pixel Pi_XOR_Pj (pixel Pi, pixel Pj)
{
    //xorarea a 2 pixeli Pi si Pj
    pixel P;
    P.B = Pi.B ^ Pj.B;
    P.G = Pi.G ^ Pj.G;
    P.R = Pi.R ^ Pj.R;
    return P;
}

pixel P_XOR_X (pixel P, unsigned int X)
{
    //xorarea a unui pixel si a unui intreg fara semn
    //X este pe 32 biti si il impatim in 4 octeti
    typedef union
    {
        unsigned int n;
        char octeti[4];
    } numar;

    numar Y;
    Y.n = X;
    pixel Pr;
    unsigned char x;
    x = Y.octeti[0];
    Pr.B = P.B ^ x;
    x = Y.octeti[1];
    Pr.G = P.G ^ x;
    x = Y.octeti[2];
    Pr.R = P.R ^ x;
    return Pr;
}

unsigned int* XORSHIFT32 (char *nume_imagine, unsigned int R0)
{
    unsigned int r, k;
    unsigned int W, H, padding;
    dimensiuni_imagine (nume_imagine, &H, &W, &padding);
    unsigned int *R = (unsigned int*) malloc (2 * H * W * sizeof (unsigned int) );
    r = R[0] = R0;
    for (k = 1; k < 2 * W * H; k++)
    {
        r = r ^ r << 13;
        r = r ^ r >> 17;
        r = r ^ r << 5;
        R[k] = r;
    }
    return R;
}

unsigned int* Algoritmul_lui_Durstenfeld (char *nume_fisier, int W, int H, int R0)
{
    //permutari aleatoare p de lungime n
    unsigned int r, k;
    unsigned int *p = (unsigned int*) malloc (H * W * sizeof (unsigned int) );
    unsigned int *R = XORSHIFT32 (nume_fisier, R0);
    for (k = 0; k < W * H; k++)
        p[k] = k;
    for (k = W * H - 1; k >= 1; k--)
    {
        r = R[k] % (k + 1);
        int aux = p[r];
        p[r] = p[k];
        p[k] = aux;
    }
    free (R);
    return p;
}

pixel* Liniarizare (char *nume_fisier)
{
    FILE *fin;
    fin = fopen (nume_fisier, "rb");
    if (fin == NULL)
    {
        printf ("Citire gresita.");
        exit (0);
    }

    //Citire header
    int i;
    unsigned char x;
    for (i = 0; i <= 53; i++)
        fread (&x, 1, 1, fin);

    unsigned int H, W, padding;
    dimensiuni_imagine (nume_fisier, &H, &W, &padding);
    pixel *L = (pixel*) malloc ( (W * H) * sizeof (pixel) );
    int j, k = 0;
    unsigned char R, G, B;
    for (i = 0; i < H; i++)
    {
        for (j = 0; j < W; j++)
        {
            fread (&B, sizeof (unsigned char), 1, fin);
            fread (&G, sizeof (unsigned char), 1, fin);
            fread (&R, sizeof (unsigned char), 1, fin);
            L[k].R = R;
            L[k].G = G;
            L[k].B = B;
            k++;
        }
        fseek (fin, padding, SEEK_CUR);
    }
    fclose (fin);
    return L;
}

char* Salvare_imagine_forma_liniarizata (char *nume_fisier)
{
    char *nume_imagine_liniarizata = (char*) malloc (30 * sizeof (char) );
    nume_imagine_liniarizata = "imagine_copiata.bmp";
    FILE *fout = NULL;
    fout = fopen (nume_imagine_liniarizata, "w+b");
    if (fout == NULL)
    {
        printf ("Citire gresita.");
        exit (0);
    }

    FILE *fin = NULL;
    fin = fopen (nume_fisier, "rb");
    if (fin == NULL)
    {
        printf ("Citire gresita.");
        exit (0);
    }

    //Citire si copiere header
    int i;
    unsigned char x;
    for (i = 0; i < 54; i++)
    {
        fread (&x, 1, 1, fin);
        fwrite (&x, 1, 1, fout);
    }
    fflush (fout);
    pixel *L = NULL;
    L = Liniarizare (nume_fisier);
    unsigned int H, W, padding;
    dimensiuni_imagine (nume_fisier, &H, &W, &padding);
    int k;
    x = 0;
    for (k = 0; k < W * H; k++)
    {
        fwrite (&L[k].B, sizeof (unsigned char), 1, fout);
        fwrite (&L[k].G, sizeof (unsigned char), 1, fout);
        fwrite (&L[k].R, sizeof (unsigned char), 1, fout);
        if (k % W == 0 && padding != 0 && k != 0)
        {
            int i;
            for (i = 1; i <= padding; i++)
            {
                fwrite (&x, sizeof (unsigned char), 1, fout);
            }
        }
        fflush (fout);
    }
    for (i = 1; i <= padding; i++)
    {
        fwrite (&x, sizeof (unsigned char), 1, fout);
    }
    fflush (fout);
    free (L);
    fclose (fin);
    fclose (fout);
    return nume_imagine_liniarizata;
}


void Permutare_pixelii_imaginii (char *nume_fisier, unsigned int H, unsigned int W, unsigned int padding, pixel **P, unsigned int R0)
{
    unsigned int *o;
    o = Algoritmul_lui_Durstenfeld (nume_fisier, W, H, R0);

    FILE *fin;
    fin = fopen (nume_fisier, "rb");
    if (fin == NULL)
    {
        printf ("Citire gresita.");
        return;
    }

    FILE *fout;
    char *nume_imagine_permutata = (char*) malloc (30 * sizeof (char) );
    nume_imagine_permutata = "imagine_permutata.bmp";
    fout = fopen (nume_imagine_permutata, "wb");
    if (fout == NULL)
    {
        printf ("Citire gresita.");
        return;
    }

    //Citire si copiere header
    int i;
    unsigned char x;
    for (i = 0; i < 54; i++)
    {
        fread (&x, 1, 1, fin);
        fwrite (&x, 1, 1, fout);
    }
    fflush (fout);
    pixel *L = Liniarizare (nume_fisier);
    x = 0;
    for (int k = 0; k < W * H; k++)
    {
        if (k % W == 0 && padding != 0 && k != 0)
        {
            int i;
            for (i = 1; i <= padding; i++)
            {
                fwrite (&x, sizeof (unsigned char), 1, fout);
            }
        }
        (*P) [k] = L[o[k]];
        fwrite (& (*P) [k].B, 1, 1, fout);
        fwrite (& (*P) [k].G, 1, 1, fout);
        fwrite (& (*P) [k].R, 1, 1, fout);
        fflush (fout);
    }
    for (i = 1; i <= padding; i++)
    {
        fwrite (&x, sizeof (unsigned char), 1, fout);
    }
    fflush (fout);
    free (o);
    free (nume_imagine_permutata);
    free (L);
    fclose (fin);
    fclose (fout);
}

pixel* Criptare (char *nume_imagine_initiala, char *nume_imagine_criptata, char *nume_fisier_text_cheie_secreta)
{
    unsigned int H, W, padding;
    dimensiuni_imagine (nume_imagine_initiala, &H, &W, &padding);
    pixel *C = (pixel*) malloc (H * W * sizeof (pixel) );
    int k;

    //Incarea valorilor R0,SV care reprezinta cheia secreta
    unsigned int SV, R0;
    FILE *in;
    in = fopen (nume_fisier_text_cheie_secreta, "r");
    if (in == NULL)
    {
        printf ("Citire gresita.");
        exit (0);
    }
    fscanf (in, "%u %u", &R0, &SV);
    fclose (in);

    unsigned int *R = XORSHIFT32 (nume_imagine_initiala, R0);
    pixel *P = (pixel*) malloc (H * W * sizeof (pixel) );
    Permutare_pixelii_imaginii (nume_imagine_initiala, H, W, padding, &P, R0);
    C[0] = P_XOR_X (P[0], SV);
    C[0] = P_XOR_X (C[0], R[W * H]);
    for (k = 1; k < W * H; k++)
    {
        C[k] = Pi_XOR_Pj (C[k - 1], P[k]);
        C[k] = P_XOR_X (C[k], R[W * H + k]);
    }

    FILE *fout;
    fout = fopen (nume_imagine_criptata, "w+b");
    if (fout == NULL)
    {
        printf ("Citire gresita.");
        exit (0);
    }

    FILE *fin;
    fin = fopen (nume_imagine_initiala, "rb");
    if (fin == NULL)
    {
        printf ("Citire gresita.");
        exit (0);
    }

    //Citire si copiere header
    int i;
    unsigned char x;
    for (i = 0; i < 54; i++)
    {
        fread (&x, 1, 1, fin);
        fwrite (&x, 1, 1, fout);
    }
    fflush (fout);

    //Scriere in imagine_criptata
    x = 0;
    for (k = 0; k < W * H; k++)
    {
        if (k % W == 0 && padding != 0 && k != 0)
        {
            int i;
            for (i = 1; i <= padding; i++)
            {
                fwrite (&x, sizeof (unsigned char), 1, fout);
            }
        }
        fwrite (&C[k].B, 1, 1, fout);
        fwrite (&C[k].G, 1, 1, fout);
        fwrite (&C[k].R, 1, 1, fout);
        fflush (fout);
    }
    for (i = 1; i <= padding; i++)
    {
        fwrite (&x, sizeof (unsigned char), 1, fout);
    }
    fflush (fout);
    free (R);
    free (P);
    fclose (fin);
    fclose (fout);
    return C;
}

void Decriptare (char *nume_imagine_initiala, char* nume_imagine_criptata, char* nume_fisier_cheie)
{
    char *nume_imagine_decriptata = (char*) malloc (30 * sizeof (char) );
    nume_imagine_decriptata = "imagine_decriptata.bmp";

    //Incarea valorilor R0,SV care reprezinta cheia secreta
    unsigned int SV, R0;
    FILE *in;
    in = fopen (nume_fisier_cheie, "r");
    if (in == NULL)
    {
        printf ("Citire gresita.");
        return;
    }
    fscanf (in, "%u %u", &R0, &SV);
    fclose (in);

    unsigned int W, H, padding;
    dimensiuni_imagine (nume_imagine_initiala, &H, &W, &padding);
    unsigned int *R = XORSHIFT32 (nume_imagine_initiala, R0);
    unsigned int *o = Algoritmul_lui_Durstenfeld (nume_imagine_initiala, W, H, R0);
    pixel *C = Liniarizare (nume_imagine_criptata);

    FILE *fin;
    fin = fopen (nume_imagine_criptata, "rb");
    if (fin == NULL)
    {
        printf ("Citire gresita.");
        return;
    }

    FILE *fout;
    fout = fopen (nume_imagine_decriptata, "wb");
    if (fout == NULL)
    {
        printf ("Citire gresita.");
        return;
    }

    //Citire si copiere header
    int i;
    unsigned char x;
    for (i = 0; i < 54; i++)
    {
        fread (&x, 1, 1, fin);
        fwrite (&x, 1, 1, fout);
    }
    fflush (fout);
    fclose (fin);

    //Inversa relatiei de substitutie folosita in procesul de criptare
    pixel *c = (pixel*) malloc (H * W * sizeof (pixel) );
    c[0] = P_XOR_X (C[0], SV);
    c[0] = P_XOR_X (c[0], R[W * H]);
    int k;
    for (k = 1; k < H * W; k++)
    {
        c[k] = Pi_XOR_Pj (C[k - 1], C[k]);
        c[k] = P_XOR_X (c[k], R[W * H + k]);
    }

    //Decriptarea imaginii
    pixel *D = (pixel*) malloc (H * W * sizeof (pixel) );
    for (k = 0; k < H * W; k++)
        D[o[k]] = c[k];

    fseek (fout, 54, SEEK_SET);
    for (k = 0; k < H * W; k++)
    {
        if (k % W == 0 && padding != 0 && k != 0)
        {
            int i;
            unsigned char x = 0;
            for (i = 1; i <= padding; i++)
            {
                fwrite (&x, sizeof (unsigned char), 1, fout);
            }
        }
        fwrite (&D[k].B, sizeof (unsigned char), 1, fout);
        fwrite (&D[k].G, sizeof (unsigned char), 1, fout);
        fwrite (&D[k].R, sizeof (unsigned char), 1, fout);
        fflush (fout);
    }
    for (i = 1; i <= padding; i++)
    {
        fwrite (&x, sizeof (unsigned char), 1, fout);
    }
    fflush (fout);
    fclose (fout);
    free (o);
    free (R);
    free (nume_imagine_decriptata);
    free (c);
    free (D);
    free (C);
}

void Testul_chi_patrat (char *nume_imagine)
{
    unsigned int *fR = (unsigned int*) malloc (256 * sizeof (unsigned int) );
    unsigned int *fG = (unsigned int*) malloc (256 * sizeof (unsigned int) );
    unsigned int *fB = (unsigned int*) malloc (256 * sizeof (unsigned int) );
    int i;
    unsigned int H, W, padding;
    dimensiuni_imagine (nume_imagine, &H, &W, &padding);
    pixel *C = Liniarizare (nume_imagine);
    for (i = 0; i < 256; i++)
    {
        fR[i] = 0;
        fG[i] = 0;
        fB[i] = 0;
    }
    for (i = 0; i < H * W; i++)
    {
        fR[C[i].R]++;
        fG[C[i].G]++;
        fB[C[i].B]++;
    }
    double F;// f cu bara deaspura
    F = (double) (H * W) / 256.0;
    double XR, XG, XB, x;
    XR = XB = XG = 0;
    for (i = 0; i < 256; i++)
    {
        x = (fR[i] - F) * (fR[i] - F) / F;
        XR = XR + x;
        x = (fG[i] - F) * (fG[i] - F) / F;
        XG = XG + x;
        x = (fB[i] - F) * (fB[i] - F) / F;
        XB = XB + x;
    }
    printf ("\n%f %f %f", XR, XG, XB);
    free (fR);
    free (fB);
    free (fG);
    free (C);
}

void grayscale_image (char* nume_fisier_sursa, char* nume_fisier_destinatie)
{
    FILE *fin, *fout;
    unsigned int dim_img;
    unsigned char pRGB[3], aux;

    fin = fopen (nume_fisier_sursa, "rb");
    if (fin == NULL)
    {
        printf ("nu am gasit imaginea sursa din care citesc");
        return;
    }

    fout = fopen (nume_fisier_destinatie, "wb+");

    fseek (fin, 2, SEEK_SET);
    fread (&dim_img, sizeof (unsigned int), 1, fin);

    //dimensiuni imagine, inaltime latime si padding
    unsigned int latime_img, inaltime_img, padding;
    dimensiuni_imagine (nume_fisier_sursa, &inaltime_img, &latime_img, &padding);

    //copiaza octet cu octet imaginea initiala in cea noua
    fseek (fin, 0, SEEK_SET);
    unsigned char c;
    while (fread (&c, 1, 1, fin) == 1)
    {
        fwrite (&c, 1, 1, fout);
        fflush (fout);
    }
    fclose (fin);

    fseek (fout, 54, SEEK_SET);
    int i, j;
    for (i = 0; i < inaltime_img; i++)
    {
        for (j = 0; j < latime_img; j++)
        {
            //citesc culorile pixelului
            fread (pRGB, 3, 1, fout);
            //fac conversia in pixel gri
            aux = 0.299 * pRGB[2] + 0.587 * pRGB[1] + 0.114 * pRGB[0];
            pRGB[0] = pRGB[1] = pRGB[2] = aux;
            fseek (fout, -3, SEEK_CUR);
            fwrite (pRGB, 3, 1, fout);
            fflush (fout);
        }
        fseek (fout, padding, SEEK_CUR);
    }
    fclose (fout);
}

//O structura in care retin x si y pt fiecare imagine care are pragul mai mare ca ps
typedef struct
{
    int Ox, Oy;
    double corelatia;
    char *template;
} coordonate;

void eliberare_matrice_pixel (pixel **matrice, unsigned int n)
{
    for (unsigned int i = 0; i < n; ++i)
        free (matrice[i]);
    free (matrice);
}

void eliberare_matrice_coordonate (coordonate **matrice, unsigned int n)
{
    for (unsigned int i = 0; i < n; ++i)
        free (matrice[i]);
    free (matrice);
}

pixel** Matrice_pixeli (char *nume_imagine)
{
    unsigned int H, W, padding;
    dimensiuni_imagine (nume_imagine, &H, &W, &padding);
    int i, j;
    pixel **M = (pixel**) malloc (H * sizeof (pixel*) );

    FILE *fin;
    fin = fopen (nume_imagine, "rb");
    if (fin == NULL)
    {
        printf ("Citire gresita.");
        return NULL;
    }
    fseek (fin, 54, SEEK_SET); //sarim de header

    //Creare matrice
    unsigned char x;
    for (i = 0; i < H; i++)
    {
        M[i] = (pixel*) malloc (W * sizeof (pixel) );
        for (j = 0; j < W; j++)
        {
            fread (&x, 1, 1, fin);
            M[i][j].B = x;
            fread (&x, 1, 1, fin);
            M[i][j].G = x;
            fread (&x, 1, 1, fin);
            M[i][j].R = x;
        }
        fseek (fin, padding, SEEK_CUR);
    }

    fclose (fin);
    return M;
}

/// O functie care face din imagine matrice, dar scrisa invers
///Avand pixelii asezati conform imaginii
pixel** Matrice_inversa (char *nume_imagine)
{
    unsigned int H, W, padding;
    dimensiuni_imagine (nume_imagine, &H, &W, &padding);
    pixel *L = Liniarizare (nume_imagine);
    pixel **Mi = (pixel**) malloc (H * sizeof (pixel*) );
    int i, j, k;
    k = H * W - 1;
    for (i = 0;  i < H;  i++)
    {
        Mi[i] = (pixel*) malloc (W * sizeof (pixel) );
        for (j = 0; j < W; j++)
        {
            Mi[i][W - j - 1].B = L[k].B;
            Mi[i][W - j - 1].G = L[k].G;
            Mi[i][W - j - 1].R = L[k].R;
            k--;
        }
    }
    return Mi;///Cred ca nu contine nimic ca altfel nu are sens
}

double corelatie (char *S, pixel **Sm, pixel **f1m, int x, int y)
{
    unsigned int H, W, padding;//Dimensiuni sablon
    dimensiuni_imagine (S, &H, &W, &padding);
    int n;//nr de pieli in sablonul S
    n = H * W;
    int i, j; //linia si coloana

    //S barat, media valorilor intensitatilor grayscale
    double Sb = 0;
    double f1b = 0;
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
        {
            Sb = Sb + Sm[i][j].R ;
            f1b = f1b + f1m[i + y][j + x].R;
        }
    Sb = Sb / n;
    f1b = f1b / n;

    //D rep derivatia standard a val intensitatilor grayscale
    double Ds, Df1;
    Ds = 0;
    Df1 = 0;
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
        {
            Ds = Ds + (Sm[i][j].R - Sb) * (Sm[i][j].R - Sb);
            Df1 = Df1 + (f1m[i + y][j + x].R - f1b) * (f1m[i + y][j + x].R - f1b) ;
        }
    Ds = sqrt (Ds / (n - 1) );
    Df1 = sqrt (Df1 / (n - 1) );

    //calculam corelatia dintre sablonul S si fereastra f1
    if (n * Ds * Df1 == 0)
        return 0;
    double corr = 0;
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            corr = corr + (f1m[i + y][j + x].R - f1b ) * (Sm[i][j].R - Sb );
    corr = corr / n / Ds / Df1;
    return corr;
}

void template_matching (char *I, char *S, double ps, int *n, coordonate **fI) //n este nr de elemente a lui fI
///fI furnizat
{
    //Crearea imaginilor alb-negru
    char *S_grayscale = (char*) malloc (30 * sizeof (char) );
    S_grayscale = "S_grayscale.bmp";
    char *I_grayscale = (char*) malloc (30 * sizeof (char) );
    I_grayscale = "I_grayscal.bmp";
    grayscale_image (I, I_grayscale);
    grayscale_image (S, S_grayscale);

    //Dimensiuni imaginii
    unsigned int H, W, padding, Hs, Ws, paddings;
    dimensiuni_imagine (S_grayscale, &Hs, &Ws, &paddings); //sablon
    dimensiuni_imagine (I_grayscale, &H, &W, &padding); //imaginea initiala


    //O structura in care retin x,y pt fiecare imagine care are pragul mai mare ca ps
    *fI = (coordonate*) malloc ( W * H * sizeof (coordonate) );

    //Calcularea corelatiei dintre sablon si fereastra(x,y)
    double C;
    pixel **M = Matrice_inversa (I_grayscale);
    pixel **Sm = Matrice_inversa (S_grayscale);
    //fereastra(x,y) din imaginea I_grayscale
    //x si y sunt coordonatele punctului de inceput al ferestrei
    int x, y;
    int k = 0;
    for (y = 0; y < H - Hs; y++)
    {
        for ( x = 0; x < W - Ws - paddings; x++)
        {
            C = corelatie (S_grayscale, Sm, M, x, y);
            if (C >= ps)
            {

                (*fI) [k].Ox = x;
                (*fI) [k].Oy = y;
                (*fI) [k].corelatia = C;
                (*fI) [k].template = (char*) malloc (30 * sizeof (char) );
                strcpy ( (*fI) [k].template, S);
                k++;
            }
        }
    }
    *n = k;
    *fI = (coordonate*) realloc (*fI, k * sizeof (coordonate) );
    eliberare_matrice_pixel (M, H);
    eliberare_matrice_pixel (Sm, Hs);
    free (I_grayscale);
    free (S_grayscale);
}

void colorare_fereastra (char *I, char *S, coordonate fereastra, pixel C)
{
    //Deschidem imaginea I pentru a putea fi editata
    FILE *fout;
    fout = fopen (I, "rb+");
    if (fout == NULL)
    {
        printf ("Citire gresita.");
        return;
    }

    int x, y;
    x = fereastra.Ox;
    y = fereastra.Oy;

    //Dimensiunile imaginii
    unsigned int H, W, padding;
    dimensiuni_imagine (I, &H, &W, &padding);
    unsigned int Hs, Ws, paddings;
    dimensiuni_imagine (S, &Hs, &Ws, &paddings);


    //Sarim de header
    fseek (fout, 54, SEEK_SET);

    //Cautam pozitia punctului de coordonate x si y in imagine
    ///Imaginea e scrisa de jos in sus

    fseek (fout, 0, SEEK_END);//coltul din dreapta sus

    //Cautam pozitia punctului de coordonate x si y in imagine
    fseek (fout, -y * 3 * W - y * padding + 3 * x - 3 * W - padding, SEEK_END);

    //Marginea superioara
    int i;
    for (i = 0; i < Ws; i++)
    {
        fwrite (&C.B, 1, 1, fout);
        fwrite (&C.G, 1, 1, fout);
        fwrite (&C.R, 1, 1, fout);
    }
    fflush (fout);
    //Marginile laterale
    //Sarim la coltul din stanga sus al ferestrei
    fseek (fout, -y * 3 * W - y * padding + 3 * x - 3 * W - padding, SEEK_END);
    for (i = 1; i <= Hs; i++)
    {
        fwrite (&C.B, 1, 1, fout);
        fwrite (&C.G, 1, 1, fout);
        fwrite (&C.R, 1, 1, fout);
        fflush (fout);
        fseek (fout, + (3 * (Ws - 1) ), SEEK_CUR);
        fwrite (&C.B, 1, 1, fout);
        fwrite (&C.G, 1, 1, fout);
        fwrite (&C.R, 1, 1, fout);
        fseek (fout, - (3 * W + padding + paddings) - 3 * Ws, SEEK_CUR);
        fflush (fout);
    }

    //Marginea superioara
    //Sarim la coltul din stanga jos al imaginii
    fseek (fout, -y * 3 * W - y * padding + 3 * x - 3 * W - padding - 3 * W * Hs - Hs * padding, SEEK_END);
    for (i = 0; i <= Ws; i++)
    {
        fwrite (&C.B, 1, 1, fout);
        fwrite (&C.G, 1, 1, fout);
        fwrite (&C.R, 1, 1, fout);
    }
    fflush (fout);
    fclose (fout);
}

coordonate* creare_vector_Detectii (int *numar_elemente)
{
    //Trebuie sa citesc numele imaginilor
    FILE *fin;
    fin = fopen ("date_program_partea2.txt", "r");
    if (fin == NULL)
    {
        printf ("Citire gresita");
        return NULL;
    }
    char *nume_imagine_initiala = (char*) malloc (30 * sizeof (char) );
    fscanf (fin, "%s", nume_imagine_initiala);
    char *sablon = (char*) malloc (30 * sizeof (char) );
    fscanf (fin, "%s", sablon);

    coordonate *fI, *D;
    int n;
    double ps = 0.5;
    template_matching (nume_imagine_initiala, sablon, ps, &n, &fI);
    D = (coordonate*) malloc (n * sizeof (coordonate) );
    int nr = 0; // Numarul de detectii
    int i;
    for (i = 0; i < n; i++)
    {
        D[nr].Ox = fI[i].Ox;
        D[nr].Oy = fI[i].Oy;
        D[nr].corelatia = fI[i].corelatia;
        D[nr].template = (char*) malloc (30 * sizeof (char) );
        strcpy (D[nr].template, fI[0].template);
        nr++;
    }
    free (fI);
    for (i = 1; i <= 9; i++)
    {
        fscanf (fin, "%s", sablon);
        template_matching (nume_imagine_initiala, sablon, ps, &n, &fI);

        D = (coordonate*) realloc (D, (nr + n) * sizeof (coordonate) );
        for (int j = 0; j < n; j++)
        {
            D[nr].Ox = fI[j].Ox;
            D[nr].Oy = fI[j].Oy;
            D[nr].corelatia = fI[j].corelatia;
            D[nr].template = (char*) malloc (30 * sizeof (char) );
            strcpy (D[nr].template, fI[i].template);
            nr++;
        }
        free (fI);
    }
    *numar_elemente = nr - 1;
    fclose (fin);
    free (nume_imagine_initiala);
    free (sablon);
    return D;
}

int cmp (const void *a, const void  *b)
{
    coordonate A, B;
    A = * (coordonate*) a;
    B = * (coordonate*) b;
    if (A.corelatia < B.corelatia)
        return 1;
    else
        return -1;
}

int Aria (int lungime, int latime)
{
    //Aria este int, deoarece lucram cu valori intregi
    return lungime * latime;
}

double suprapunere (coordonate di, coordonate dj,unsigned int Hi, unsigned int Wi)
{
    //di si dj sunt coordonatele punctelor de inceput(colt stanga sus) a ferestrelor i si j
    double Arie_di, Arie_dj;
    Arie_di=Arie_dj= Aria (Hi, Wi);
    double Arie_intersectie;
    if (abs (dj.Oy - di.Oy) <= Hi && abs (dj.Ox - di.Ox) <= Wi)
    {
        int lungime, latime;
        latime = Wi - abs (di.Ox - dj.Ox);
        lungime = Hi - abs (di.Oy - dj.Oy);
        Arie_intersectie = Aria (lungime, latime);
    }
    else
        Arie_intersectie = 0;

    double aux = Arie_di + Arie_dj - Arie_intersectie;
    return Arie_intersectie / aux;
}

coordonate* eliminare_non_maxime (coordonate *D, int nr, int *nr_final)
{
    int i, j;
    double supra;
    coordonate *F = (coordonate*) malloc (nr * sizeof (coordonate) );
    int k = 0;
    int ok = 0;

    unsigned int Hs, Ws, paddings;
    dimensiuni_imagine (D[1].template, &Hs, &Ws, &paddings);

    for (i = nr - 1; i > 0; i--)
    {
        ok = 1;
        for (j = 0; j < i - 1; j++)
        {
            //Corelatie D[i]<D[j] si i>j
            supra = suprapunere (D[i], D[j], Hs, Ws);
            if (supra >= 0.2) //Ferestrele se suprapun
            {
                ok = 0;
            }
        }
        if (ok == 1)
        {
            F[k].Ox = D[i].Ox;
            F[k].Oy = D[i].Oy;
            F[k].corelatia = D[i].corelatia;
            F[k].template = (char*) malloc (30 * sizeof (char) );
            strcpy (F[k].template, D[i].template);
            k++;
        }
    }
    *nr_final = k;
    free (D);
    F = (coordonate*) realloc (F, k * sizeof (coordonate) );
    return F;
}

int main()
{

    //Prima parte
    FILE *fin;
    fin = fopen ("date_program_partea1.txt", "r");
    if (fin == NULL)
    {
        printf ("Nu se poate deschide fisierul text");
        return 0;
    }
    //Criptare si Decriptare
    char *nume_imagine_initiala = (char*) malloc (30 * sizeof (char) );
    char *nume_imagine_criptata = (char*) malloc (30 * sizeof (char) );
    char *nume_fisier_text_cheie_secreta = (char*) malloc (30 * sizeof (char) );
    fscanf (fin, "%s %s %s", nume_imagine_initiala, nume_imagine_criptata, nume_fisier_text_cheie_secreta);
    Criptare (nume_imagine_initiala, nume_imagine_criptata, nume_fisier_text_cheie_secreta);
    Decriptare (nume_imagine_initiala, nume_imagine_criptata, nume_fisier_text_cheie_secreta);
    printf ("Valorile testului chi patrat pentru imaginea initiala: ");
    Testul_chi_patrat (nume_imagine_initiala);
    printf ("\nValorile testului chi patrat pentru imaginea criptata: ");
    Testul_chi_patrat (nume_imagine_criptata);
    printf ("\n");
    fclose (fin);

    //A 2-a parte
    FILE *in;
    in = fopen ("date_program_partea2.txt", "r");
    if (fin == NULL)
    {
        printf ("Nu se poate deschide fisierul text");
        return 0;
    }
    char *nume_imagine = (char*) malloc (30 * sizeof (char) );
    fscanf (in, "%s", nume_imagine);
    char **sablon = (char**) malloc (10 * sizeof (char*) );
    int j;
    for (j = 0; j <= 9; j++)
    {
        sablon[j] = (char*) malloc (30 * sizeof (char) );
        fscanf (in, "%s", sablon[j]);
    }
    int nr, nr_elemente;
    coordonate *D = creare_vector_Detectii (&nr);
    qsort (D, nr, sizeof (coordonate), cmp);
    coordonate *F = eliminare_non_maxime (D, nr, &nr_elemente);

    //Culori pentru bordare
    pixel *C = (pixel*) malloc (10 * sizeof (pixel) );
    C[0].R = 255;
    C[0].G = C[0].B = 0;
    C[1].R = C[1].G = 255;
    C[1].B = 0;
    C[2].R = C[2].B = 0;
    C[2].G = 255;
    C[3].R = 0;
    C[3].G = C[3].B = 255;
    C[4].R = C[4].B = 255;
    C[4].G = 0;
    C[5].R = C[5].G = 0;
    C[5].B = 255;
    C[6].R = C[6].G = C[6].B = 192;
    C[7].R = 255;
    C[7].G = 140;
    C[7].B = 0;
    C[8].R = C[8].B = 128;
    C[8].G = 0;
    C[9].R = 128;
    C[9].G = C[9].B = 0;

    ///Copiez imaginea test.bmp si calorez intr o copie
    char *nume_imagine_colorata = (char*) malloc (30 * sizeof (char) );
    nume_imagine_colorata = Salvare_imagine_forma_liniarizata (nume_imagine);
    int i;
    for (i = nr_elemente - 1; i >= 0; i--)
    {
        for (j = 0; j <= 9; j++)
            if (strcmp (F[i].template, sablon[j]) == 0)
                break;
        colorare_fereastra (nume_imagine_colorata, F[i].template, F[i], C[j]);
    }
    fclose (in);

    //Eliberare memorie alocata dinamic
    free (D);
    free (F);
    free (nume_fisier_text_cheie_secreta);
    free (nume_imagine);
    free (nume_imagine_colorata);
    free (nume_imagine_criptata);
    free (nume_imagine_initiala);
    free (C);
    return 0;
}
