unsigned short* read(std::string filename, std::map<std::string, std::string>& pars)
{
    std::ifstream fin(filename, std::ios::binary);
    if (!fin.is_open())
        return NULL;

    std::string name, value;
    for (int i = 0; i < 8; i++)
    {
    	fin >> name >> value;
	std::cout << name << " " << value << std::endl;
	pars[name] = value;
    }
    char ch;
    fin.get(ch);

    int cols = std::stoi(pars["xSize"]);
    int rows = std::stoi(pars["ySize"]);

    unsigned short* image = new unsigned short[cols * rows];
    fin.read((char*)(image), sizeof(image[0])*(cols * rows));
    //for (int i = 0; i < cols * rows; i++)
    //    fin >> image[i];

    std::cout << image[100] << " " << image[667] << std::endl;
    return image;
}

TH2I* read_and_fill(std::string filename, double factor = 10)
{
    std::map<std::string, std::string> pars;
    unsigned short* image = read(filename, pars);
    int cols = std::stoi(pars["xSize"]);
    int rows = std::stoi(pars["ySize"]);
    static int num = 0;
    /*
    TH2I* hist2d = new TH2I(Form("hist%d", num), Form("hist%d", num++), cols, 0, cols, rows, 0, rows);
    for (int j = 0; j < rows; j++)
        for (int i = 0; i < cols; i++)
            hist2d->SetBinContent(i, rows - j - 1, image[(j * cols) + i]);
    */

    TH2I* hist2d = new TH2I(Form("hist%d", num), Form("hist%d", num++), rows, 0, 3.69 * rows * 1e-3 * factor, cols, 0, 3.69 * cols * 1e-3 * factor);
    for (int j = 0; j < rows; j++)
        for (int i = 0; i < cols; i++)
            hist2d->SetBinContent(j, i, image[(j * cols) + i]);

    hist2d->GetXaxis()->SetTitle("<- Beam direction <-, mm");
    hist2d->GetYaxis()->SetTitle("Vertical, mm");

    delete[] image;

    return hist2d;
}

TH1I* build_1dimhist(TH2I* hist2d)
{
    //TFile* f = TFile::Open(filename);
    //TH2I* hist2d = (TH2I*)f->Get("hist2d");
    int* data = hist2d->GetArray();

    TH1I* hist1d = new TH1I("hist1d", "hist1d", 65536, 0, 65536);

    for (int i = 0; i < hist2d->GetEntries(); i++)
        hist1d->Fill(data[i]);

    //hist1d->Draw();
    return hist1d;
}

std::pair<std::vector<double>, std::vector<double>> calc_stat(const std::vector<std::vector<int>> array)
{
	std::vector<double> means;
	means.resize(3388*2712, 0);
	std::vector<double> devs;
	devs.resize(3388*2712, 0);

	int vsize = array[0].size();
	for (int i = 0; i < vsize; i++)
	{
	    int count = 0;
	    for (int j = 0; j < array.size(); j++)
	    {
		if (array.at(j).at(i) > 40000)
		    continue;

	        means.at(i) += array.at(j).at(i);
		count++;
	    }

	    means.at(i) *= 1. / count;
	}

	for (int i = 0; i < vsize; i++)
	{
	    int count = 0;
            for (int j = 0; j < array.size(); j++)
	    {
		if (array.at(j).at(i) > 40000)
		    continue;

	        devs.at(i) += pow(means.at(i) - array.at(j).at(i), 2);
		count++;
	    }

	    devs.at(i) *= 1./count;
	}

	for (int i = 0; i < 3388*2712; i++)
	    devs.at(i) = sqrt(devs.at(i));

	return {means, devs};
}

std::vector<int> fill_data(TString fname, double& time, bool subtract_ped = true)
{
    //TFile* f = TFile::Open(fname);
    //TH2I* hist2d = (TH2I*)f->Get("hist2d");
    //TParameter<double>* expTime = (TParameter<double>*)f->Get("exposureTime");
    //time = expTime->GetVal();

    std::map<std::string, std::string> pars;
    unsigned short* ped = read("data/photo_2025-01-22T02:15:39.042.dat", pars);
    unsigned short* image = read(std::string(fname.Data()), pars);

    int cols = std::stoi(pars["xSize"]);
    int rows = std::stoi(pars["ySize"]);
    // Image size 3388 x 2712
    //int cols = 3388;
    //int rows = 2712;
    std::vector<int> res;
    res.reserve(cols * rows);

    /*for (int i = 0; i < cols; i++){
        for (int j = 0; j < rows; j++)
	{
            res.push_back(hist2d->GetBinContent(i, j));
	}
    }*/

    for (int i = 0; i < cols * rows; i++)
        res.push_back(image[i] - ped[i]);

    time = std::stod(pars["exposureTime"]);

    return res;
}

void fill_tree(const char* dirname, const char* ext, bool subtract_ped = true)
{
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
	std::map<int, std::vector<std::vector<int>>> array;
	double time;

        while ((file=(TSystemFile*)next())) {
            fname = TString(dirname) + "/" + file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
	            std::cout << fname << std::endl;
		    auto data = fill_data(fname, time, subtract_ped);
		    array[time*1e3].push_back(data);
            }
        }

	TFile *res_file = new TFile("res.root", "UPDATE");
	TTree* tree = new TTree("tree", "tree");
	int col, row;
	double mean, dev;
	tree->Branch("col", &col);
	tree->Branch("row", &row);
	tree->Branch("mean", &mean);
	tree->Branch("dev", &dev);
	tree->Branch("time", &time);

	TH2I* hist2d = new TH2I("hist2d", "hist2d", 3388, 0, 3388, 2712, 0, 2712);
	TGraph* profileX = new TGraph();
	TGraph* profileY = new TGraph();
	profileX->SetName("profileX");
	profileY->SetName("profileY");

	for (auto it = array.begin(); it != array.end(); ++it)
	{
	    time = it->first * 1.e-3;
	    auto [means, devs] = calc_stat(it->second);
	    std::cout << means[0] << " " << devs[0] << " " << time << std::endl;
	    for (int i = 0; i < 3388; i++)
	    {
	        for (int j = 0; j < 2712; j++)
	        {
		        row = j;
		        col = i;
		        mean = means.at(3388*j + i);
		        dev = devs.at(3388*j + i);
			hist2d->SetBinContent(i, j, mean);
		        tree->Fill();
	        }
	    }
	}

	for (int ix = 0; ix < 3388; ix++)
	{
	    double ymean = 0;
	    for (int iy = 0; iy < 2712; iy++)
	    {
	        double c = hist2d->GetBinContent(ix, iy);
		ymean += c * 1./2712;
	    }
	    profileX->AddPoint(ix, ymean);
	}

	for (int iy = 0; iy < 2712; iy++)
	{
	    double xmean = 0;
	    for (int ix = 0; ix < 3388; ix++)
	    {
	        double c = hist2d->GetBinContent(ix, iy);
		xmean += c * 1./3388;
	    }
	    profileY->AddPoint(iy, xmean);
	}

	tree->Write();
	hist2d->Write();
	profileX->Write();
	profileY->Write();
    }
}

int* average(int* input, int xsize = 2, int ysize = 2)
{
    int cols = 2712;
    int rows = 3388;
    int* avg = new int[cols * rows];
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
	{
	    avg[i * cols + j] = 0;

	    if (i < ysize || i >= rows - ysize)
	        continue;

	    if (j < xsize || j >= cols - xsize)
	        continue;

	    for (int k = -ysize; k < ysize+1; k++)
	        for (int l = -xsize; l < xsize+1; l++)
		    avg[i * cols + j] += input[(i + k) * cols + (j + l)];
	    
            avg[i * cols + j] /= (2*xsize+1)*(2*ysize+1);
	}
    }

    return avg;
}

int* filter(int* input, int xsize = 4, int ysize = 4)
{
    int cols = 2712;
    int rows = 3388;
    int* fil = new int[cols * rows];
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
	    fil[i * cols + j] = input[i * cols + j] < 0 ? -1 : 1;

    int* fil2 = new int[cols * rows];
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
	{
	    if (i < ysize || i >= rows - ysize)
	    {
	        fil2[i * cols + j] = -1;
	        continue;
	    }

	    if (j < xsize || j >= cols - xsize)
	    {
	        fil2[i * cols + j] = -1;
	        continue;
	    }

	    double sum = 0;
	    for (int k = -ysize; k < ysize+1; k++)
	        for (int l = -xsize; l < xsize+1; l++)
		    sum += fil[(i + k) * cols + (j + l)];
	    
	    sum /= (2*xsize+1)*(2*ysize+1);
	    if (abs(sum)<0.7)
	        fil2[i * cols + j] = -1;
	    else
	        fil2[i * cols + j] = 1;
	}
    }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
	{
	    if (i < ysize || i >= rows - ysize)
	    {
	        fil[i * cols + j] = -1;
	        continue;
	    }

	    if (j < xsize || j >= cols - xsize)
	    {
	        fil[i * cols + j] = -1;
	        continue;
	    }

	    double sum = 0;
	    for (int k = -ysize; k < ysize+1; k++)
	        for (int l = -xsize; l < xsize+1; l++)
		    sum += fil2[(i + k) * cols + (j + l)];

	    if (sum < 0)
	        fil[i * cols + j] = -1;
	    else
	        fil[i * cols + j] = 1;
	}
    }

    delete[] fil2;
    return fil;
}


int* ptr;
void fcn(int &npar, double* gin, double& f, double* par, int iflag)
{
    f = 0;
    double x0 = par[0];
    double y0 = par[1];
    double a = par[2];
    double b = par[3];
    int cols = 2712;
    int rows = 3388;
    int count = 0;
    for (int irow = 0; irow < rows; irow++)
    {
        for (int jcol = 0; jcol < cols; jcol++)
        {
            if (ptr[irow * cols + jcol] != 1)
	        continue;

	    double x = 3.69 * jcol;
	    double y = 3.69 * irow;
	    double dist = sqrt(pow((x - x0)/a, 2) + pow((y - y0)/b, 2));
	    f += abs(dist - 1);
	    count++;
        }
    }
    f /= count;
}

double DrawEdges(std::string filename)
{
    using namespace std::placeholders;

    std::map<std::string, std::string> pars;
    unsigned short* in = read(filename, pars);

    int cols = 2712;
    int rows = 3388;

    int* transpose = new int[cols * rows];
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
	    transpose[i * cols + j] = in[j * rows + i];

    int* avg = average(transpose);
    delete[] transpose;

    int* outx = new int[cols * rows];
    int* outy = new int[cols * rows];
    double* out = new double[cols * rows];
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
	{
	    outx[i * cols + j] = 0;
	    outy[i * cols + j] = 0;
	    out[i * cols + j] = 0;
	}
    }

    for (int irow = 1; irow < rows - 1; irow++)
    {
        for (int jcol = 1; jcol < cols - 1; jcol++)
	{
	    for (int k = -1; k < 2; k++)
	    {
	        for(int l = -1; l < 2; l++)
		{
		    int gx = l;
		    int gy = k;
		    outx[irow * cols + jcol] += avg[(irow + k) * cols + (jcol + l)] * gx;
		    outy[irow * cols + jcol] += avg[(irow + k) * cols + (jcol + l)] * gy;
		}
	    }

	    /*
	    out[irow * cols + jcol] = sqrt(pow(1.*outx[irow * cols + jcol], 2) + pow(1.*outy[irow * cols + jcol], 2));
	    if (out[irow * cols + jcol] < threshold)
	        out[irow * cols + jcol] = 1;
	    */
	}
    }
    delete[] avg;

    int* avgx = average(outx, 1, 1);
    int* filtx = filter(avgx);

    TCanvas* c1 = new TCanvas("cEllipse", "cEllipse", 900, 900);
    TH2I* hist1 = new TH2I("histEllipse", "XEdges; microns; microns", cols, 0, 3.69*cols, rows, 0, 3.69*rows);
    for (int j = 0; j < rows; j++)
        for (int i = 0; i < cols; i++)
            hist1->SetBinContent(i, j, filtx[(j * cols) + i]);
    hist1->SetMaximum(1);
    hist1->SetMinimum(-1);
    hist1->Draw("colz");

    delete[] outx;
    delete[] outy;
    delete[] out;
    delete[] avgx;

    //fit
    TMinuit minuit(4);
    ptr = filtx;
    minuit.SetFCN(fcn);
    minuit.DefineParameter(0, "x0", 3.69*cols*0.5, 10, 0, 3.69*cols);
    minuit.DefineParameter(1, "y0", 3.69*rows*0.5, 10, 0, 3.69*rows);
    minuit.DefineParameter(2, "a", 3700, 10, 3000, 6000);
    minuit.DefineParameter(3, "b", sqrt(2.)*3700, 10, 3000, 9000);
    minuit.Migrad();
    double p[4];
    double ep[4];
    for (int i = 0; i < 4; i++)
    {
        minuit.GetParameter(i, p[i], ep[i]);
        std::cout << p[i] << " +- " << ep[i] << std::endl;
    }

    TGraph* gr = new TGraph();
    for (int i = 0; i < 1000; i ++)
    {
	double x = p[0] + p[2] * cos(i/1000. * 2*M_PI);
	double y = p[1] + p[3] * sin(i/1000. * 2*M_PI);
        gr->AddPoint(x, y);
    }
    gr->SetLineColor(kRed);
    gr->Draw("L""same");

    delete[] filtx;
    return 120./(p[3] * 2 * 1e-3);
}

double fline(double* x, double* par)
{
    if (x[0] > 10 && x[0] < 96)
    {
        TF1::RejectPoint();
	return 0;
    }
    return par[0] + par[1] * x[0];
}

void DrawCanvas(std::string filename, double factor = 10.5)
{
    gStyle->SetPalette(kGreyScale);
    gStyle->SetOptStat(11);
    gStyle->SetOptFit(111);

    TH2I* h = read_and_fill(filename, factor);
    h->SetTitle("Photo");
    h->SetMaximum(4000);

    /*
    int cols = 3388, rows = 2712;
    TFile* f = TFile::Open("dataset1/photo_88.root");
    TH2I* hped = (TH2I*)f->Get("hist2d");
    TH2I* hp = new TH2I("hp", "hp", cols, 0, cols, rows, 0, rows);

    for (int j = 0; j < rows; j++)
        for (int i = 0; i < cols; i++)
            hp->SetBinContent(i, rows - j - 1, hped->GetBinContent(i, j));

    //hp->SetMaximum(4000);
    //hp->Draw("colz");

    for (int j = 0; j < rows; j++)
        for (int i = 0; i < cols; i++)
	    h->SetBinContent(i, j, h->GetBinContent(i, j) - hp->GetBinContent(i, j));
    */
	    
    h->SetMinimum(0);
    h->SetMaximum(5000);


    TH1D* px = h->ProjectionX("px", 1800, 1900);
    px->SetTitle(Form("Projection X, Y = [%.2f, %.2f] mm", 3.69 * 1.8 * factor, 3.69 * 1.9 * factor));
    TH1D* py = h->ProjectionY("py", 1300, 1400);
    py->SetTitle(Form("Projection Y, X = [%.2f, %.2f] mm", 3.69 * 1.3 * factor, 3.69 * 1.4 * factor));

    TCanvas* c = new TCanvas("cPorfile", "cProfile", 1200, 900);
    TPad* p1 = new TPad("p1", "p1", 0., 0., 0.7, 1.);
    p1->Draw();
    TPad* p2 = new TPad("p2", "p2", 0.7, 0.5, 1., 1.);
    p2->Draw();
    TPad* p3 = new TPad("p3", "p3", 0.7, 0., 1., 0.5);
    p3->Draw();

    p1->cd();
    h->Draw("colz");
    p2->cd();
    px->Draw();
    p2->cd()->SetGrid();

    TF1* fxb = new TF1("fxb", fline, 0, 104, 2);
    fxb->SetParameter(1e5, 2*1e3);
    px->Fit(fxb, "RM", "", 0, 104);

    TF1* fx = new TF1("fx", "[0] + [1]*x + gaus(2)", 20, 87);
    fx->SetParNames("b", "k", "A", "mu", "sigma");
    fx->SetParameters(fxb->GetParameter(0), fxb->GetParameter(1), 111, 50, 30);
    fx->FixParameter(0, fxb->GetParameter(0));//fx->SetParLimits(0, 1e4, 2e6);
    fx->FixParameter(1, fxb->GetParameter(1));//fx->SetParLimits(1, 1e3, 2e3);
    fx->SetParLimits(2, 0, 2e9);
    fx->SetParLimits(3, 30, 90);
    fx->SetParLimits(4, 0, 100);
    px->Fit(fx, "RM", "", 20, 87);

    TF1* fx2 = new TF1("fx2", "[0] + [1]*x", 0, 104);
    fx2->SetParameters(fxb->GetParameter(0), fxb->GetParameter(1));
    fx2->Draw("same");

    p3->cd();
    py->Draw();
    p3->cd()->SetGrid();

    TF1* fy = new TF1("fy", "[0] + [1]*x + gaus(2)", 0, 125);
    fy->SetParNames("b", "k", "A", "mu", "sigma");
    fy->SetParameters(100, 110, 111, 70, 111);
    fy->SetParLimits(0, 0, 1e9);
    fy->SetParLimits(1, 0, 1e9);
    fy->SetParLimits(2, 0, 1e9);
    fy->SetParLimits(4, 0, 100);
    py->Fit(fy, "RM", "", 0, 125);

    TF1* fy2 = new TF1("fy2", "[0] + [1]*x", 0, 125);
    fy2->SetParameters(fy->GetParameter(0), fy->GetParameter(1));
    fy2->Draw("same");
}
