void Color_def(){//Model Color Palette for Color Universal Design GUIDE BOOK
	// Accent Color
	gROOT->GetColor( 1)->SetRGB( 64/255., 64/255., 64/255.);//Dark Grey
	gROOT->GetColor( 2)->SetRGB(255/255., 75/255.,  0/255.);//Red
	gROOT->GetColor( 3)->SetRGB(  0/255., 90/255.,255/255.);//Blue
	gROOT->GetColor( 4)->SetRGB(  3/255.,175/255.,122/255.);//Green
	gROOT->GetColor( 5)->SetRGB(246/255.,170/255.,  0/255.);//Orange
	gROOT->GetColor( 6)->SetRGB( 77/255.,196/255.,255/255.);//Sky
	gROOT->GetColor( 7)->SetRGB(153/255.,  0/255.,153/255.);//Purple
	gROOT->GetColor( 8)->SetRGB(255/255.,128/255.,130/255.);//Pink
	gROOT->GetColor( 9)->SetRGB(128/255., 64/255.,  0/255.);//Brown
	gROOT->GetColor(10)->SetRGB(255/255.,241/255.,  0/255.);//Yellow
	// Base Color
	gROOT->GetColor(11)->SetRGB(255/255.,202/255.,191/255.);//Light Pink
	gROOT->GetColor(12)->SetRGB(255/255.,255/255.,128/255.);//Cream
	gROOT->GetColor(13)->SetRGB(216/255.,242/255., 85/255.);//Light Yellow Green
	gROOT->GetColor(14)->SetRGB(191/255.,228/255.,255/255.);//Light Sky
	gROOT->GetColor(15)->SetRGB(255/255.,202/255.,128/255.);//Beige
	gROOT->GetColor(16)->SetRGB(119/255.,217/255.,168/255.);//Light Green
	gROOT->GetColor(17)->SetRGB(201/255.,172/255.,230/255.);//Light Purple
	// Grey
	gROOT->GetColor(18)->SetRGB(200/255.,200/255.,203/255.);//Light Grey
	gROOT->GetColor(19)->SetRGB(132/255.,145/255.,158/255.);//Grey
	gROOT->GetColor(20)->SetRGB(  0/255.,  0/255.,  0/255.);//Black
}
void Color(){
	Color_def();
	TCanvas *c = new TCanvas("c","colors",0,0,800,600);
	c->DrawColorTable();
}
