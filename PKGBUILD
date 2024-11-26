# Maintainer: Carlos Sanchez <research@sanchezcarlosjr.com>

pkgname=slc-bin
pkgver=2.3.0
pkgrel=1
pkgdesc="Precompiled TAOS II event simulator for astronomical occultation events"
arch=('x86_64')
url="https://github.com/sanchezcarlosjr/occultation_light_curves"
license=('MIT')
depends=('fftw' 'gsl' 'hdf5')
source_x86_64=("https://github.com/sanchezcarlosjr/occultation_light_curves/releases/download/v${pkgver}/slc")
sha256sums_x86_64=('7b3ef89fdc5f6cc27e6d974efb186577a951d1c39a2d0a236cc11b2ab679f995')

package() {
  install -Dm755 "slc" "${pkgdir}/usr/bin/slc"
}
