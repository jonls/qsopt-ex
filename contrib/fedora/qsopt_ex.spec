Summary: An exact, linear programming solver
Name: qsopt_ex
Version: 2.5.10_p1
Release: 1%{dist}
License: GPLv3+
Group: Applications/Engineering
URL: https://github.com/jonls/qsopt-ex
Source0: https://github.com/jonls/qsopt-ex/releases/download/v%{version}/%{name}-%{version}.tar.xz
BuildRequires: chrpath
BuildRequires: libtool
BuildRequires: zlib-devel
BuildRequires: bzip2-devel
BuildRequires: gmp-devel

%description
Exact linear programming solver. This is a fork of QSopt_ex by Daniel
Espinoza et al. version 2.5.10. The goal of this fork is to update the
software, and in particular the build system, to be more friendly. In
addition the external dependencies have been reduced by removing the
dependency on EGlib and GNU awk. The dependencies may be further reduced
later.

This is the base library.

%package -n %{name}-devel
Summary: Development files for qsopt_ex
Group: Applications/Engineering
Requires: %{name} = %{version}-%{release}

%description -n %{name}-devel
Development files for qsopt_ex, an exact linear programming solver.

%package -n %{name}-python
Summary: Python module for qsopt_ex
Group: Applications/Engineering
BuildRequires: python2-devel >= 2.7
BuildRequires: Cython >= 0.20
Requires: %{name} = %{version}-%{release}

%description -n %{name}-python
Python module for qsopt_ex, an exact linear programming solver.

%prep
%setup -q

%build
%configure --disable-static --enable-python-module
make %{?_smp_mflags} V=1

%install
rm -rf %{buildroot}
make DESTDIR=%{buildroot} install INSTALL="install -p"
chrpath --delete %{buildroot}%{_bindir}/esolver
chrpath --delete %{buildroot}%{python2_sitearch}/qsoptex.so
rm -f %{buildroot}%{_libdir}/libqsopt_ex.la
rm -f %{buildroot}%{python2_sitearch}/qsoptex.la

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%doc README.md NEWS.md License.txt
%{_bindir}/esolver
%{_libdir}/libqsopt_ex.so.0
%{_libdir}/libqsopt_ex.so.0.1.1

%files -n %{name}-devel
%{_libdir}/libqsopt_ex.so
%{_includedir}/qsopt_ex/

%files -n %{name}-python
%{python2_sitearch}/qsoptex.so

%changelog
* Sat Nov 29 2014 Jon Lund Steffensen <jonlst@gmail.com> - 
- Initial build.
