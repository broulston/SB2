import importlib

import makeSB2temp

importlib.reload(makeSB2temp)

SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M1",
                      "Lspect050g78n.dat", "WD4")
SB2.plotCompositeSB2(saveplot=False)

SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M1",
                      "Lspect050g78n.dat", "WD4")
SB2.plotCompositeSB2(saveplot=False)

SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M1",
                      "Lspect050g78n.dat", "WD4")
SB2.plotCompositeSB2(saveplot=False)

SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M1",
                      "Lspect050g78n.dat", "WD4")
SB2.plotCompositeSB2(saveplot=False)

SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M1",
                      "Lspect050g78n.dat", "WD4")
SB2.plotCompositeSB2(saveplot=False)

SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M1",
                      "Lspect050g78n.dat", "WD4")
SB2.plotCompositeSB2(saveplot=False)


def test_makeSB2():
    SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M1",
                          "143.03409+026.79849_57474_4-2071.txt", "WD5")
    SB2.plotCompositeSB2()
    SB2.saveCompositeSB2()


SB2 = makeSB2temp.SB2("142.58623+048.87419_57851_3-149693286.txt", "M4",
                      "177.14002+001.48314_57358_27-1853.txt", "B6")
SB2.plotCompositeSB2()

SB2 = makeSB2temp.SB2("159.21672+025.58140_57852_3-125612648.txt", "M0",
                      "191.10181+012.82911_57059_4-20296.txt", "K7")
SB2.plotCompositeSB2()

SB2 = makeSB2temp.SB2("142.58623+048.87419_57851_3-149693286.txt", "M4",
                      "177.14002+001.48314_57002_27-1853.txt", "B6")
SB2.plotCompositeSB2(filename="M4+B6_2.eps")

SB2 = makeSB2temp.SB2("142.58623+048.87419_57851_3-149693286.txt", "M4",
                      "177.14002+001.48314_57390_27-1853.txt", "B6")
SB2.plotCompositeSB2(filename="M4+B6_3.eps")


dist2 = 362 * u.pc

back_flux = SB2.flux / (4.0*np.pi*SB2.wavelength*dist2.to(u.cm)**2)

plt.plot(SB2.wavelength, back_flux)
plt.show()
plt.clf()
plt.close()

composite_dir = "data/composites/"
composite_plotdir = "data/composites/plotsL4/"

# *** dM + WD ***
SB2 = makeSB2temp.SB2("139.46845+021.74163_57738_7-5464872.txt", "M0",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"M0+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"M0+WD.txt")

SB2 = makeSB2temp.SB2("227.12491+036.19676_57833_3-148660775.txt", "M1",
                     "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"M1+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"M1+WD.txt")

SB2 = makeSB2temp.SB2("141.78483+055.58795_57712_3-152400105.txt", "M2",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"M2+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"M2+WD.txt")

SB2 = makeSB2temp.SB2("226.54240+035.73449_57833_3-148270084.txt", "M3",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"M3+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"M3+WD.txt")

SB2 = makeSB2temp.SB2("226.76735+035.81285_57860_3-148400110.txt", "M4",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"M4+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"M4+WD.txt")

# *** K + WD ***
SB2 = makeSB2temp.SB2("137.54807+021.89058_57058_4-4333.txt", "K2",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"K2+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"K2+WD.txt")

SB2 = makeSB2temp.SB2("190.99932-003.46924_57447_3-109319870.txt", "K3",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"K3+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"K3+WD.txt")

SB2 = makeSB2temp.SB2("161.98968-002.75814_57435_3-110441670.txt", "K4",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"K4+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"K4+WD.txt")

SB2 = makeSB2temp.SB2("176.62262+001.35084_57058_4-19883.txt", "K5",
                      "Lspect040g78n.dat", "WD3")
SB2.plotCompositeSB2(filename=composite_plotdir+"K5+WD3.eps")
#SB2.saveCompositeSB2(filename=composite_dir+"K5+WD.txt")

