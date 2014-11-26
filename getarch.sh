#!/bin/bash
# OpenEye Scientific Software
# Copyright (C) 2007, 2008, 2009, 2011, 2012, 2013
#

ostype="unknown"
osver=`uname -r | sed -n -e 's/\([0-9]\{1,2\}\.[0-9]\).*/\\1/p'`
osmaj=`echo $osver | cut -d. -f1`
osmin=`echo $osver | cut -d. -f2`
osrep=""
dist=""
arch=""

case `uname -s` in
  CYGWIN*)
    ostype="microsoft"
    case $PROCESSOR_ARCHITEW6432 in
      *64) osver="win64"
           arch="x64"
           ;;
      *)   osver="win32"
           arch="x86"
           ;;
    esac
    ;;

  Darwin)
    ostype="osx"
    case $osmaj in
      7)  osver="10.3" ;;
      8)  osver="10.4" ;;
      9)  osver="10.5" ;;
     10)  osver="10.6"; arch="x64" ;;
     11)  osver="10.7"; arch="x64" ;;
     12)  osver="10.8"; arch="x64" ;;
     13)  osver="10.9"; arch="x64" ;;
    esac

    if test "$arch" = "" ; then
      case `machine` in
        pentium4|*86) arch="x86"   ;;
        *)            arch="ppc32" ;;
      esac
    fi

    osrep="$ostype-$osver-$arch"
    ;;

  Linux)
    ostype="linux"
    osmaj=""
    osmin=""

    if test -r /etc/sgi-release ; then
      version=`grep Linux /etc/sgi-release`
      if test "$version" = "SGI ProPack 4 for Linux, Build 400r1-0502032344" ; then
        ostype="sgipropack"
        osver="4.0"
      elif test "$version" = "SGI ProPack 3SP2 for Linux, Build 302rp04110909_10116-0411090951" ; then
        ostype="sgipropack"
        osver="3.2"
      elif test "$version" = "SGI ProPack 2.3 for Linux, Build 230rp03111013_10029-0311101317" ; then
        ostype="sgipropack"
        osver="2.3"
      fi
    elif test -r /etc/redhat-release ; then
      ostype="redhat"
      version=`grep release /etc/redhat-release`
      if test "$version" = "Red Hat Linux release 9 (Shrike)" ; then
        osver="9.0"
      elif test "$version" = "Red Hat Linux release 8.0 (Psyche)" ; then
        osver="8.0"
      elif test "$version" = "Red Hat Linux release 7.3 (Valhalla)" ; then
        osver="7.3"
      elif test "$version" = "Red Hat Linux release 7.2 (Enigma)" ; then
        osver="7.2"
      elif test "$version" = "Red Hat Linux release 7.1 (Seawolf)" ; then
        osver="7.1"
      elif test "$version" = "Red Hat Linux release 7.0 (Guinness)" ; then
        osver="7.0"
      elif test "$version" = "Red Hat Linux release 6.2 (Zoot)" ; then
        osver="6.2"
      elif test "$version" = "Red Hat Linux release 6.1 (Cartman)" ; then
        osver="6.1"
      elif test "$version" = "Red Hat Linux release 6.0 (Hedwig)" ; then
        osver="6.0"
      else
        osver=`echo "$version" | sed "s/.*release\ \([^\ \.]*\)[\ \.].*/\\1/"`
        fedora=`echo "$version" | grep Fedora`
        if test "$fedora" ; then
          osver="Fedora$osver"
        else
          osver="RHEL$osver"
        fi
      fi
    elif test -r /etc/SuSE-release ; then
      ostype="suse"
      version=`grep -i suse /etc/SuSE-release`
      suse_ed=`echo "$version" | grep -i "Enterprise Desktop"`
      suse_es=`echo "$version" | grep -i "Enterprise Server"`
      suse_sles=`echo "$version" | grep -i "SLES"`
      suse_open=`echo "$version" | grep -i "open"`
      if test "$suse_ed" ; then
        osver=`echo "$version" | sed "s/.*Desktop\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      elif test "$suse_es" ; then
        osver=`echo "$version" | sed "s/.*Server\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      elif test "$suse_sles" ; then
        osver=`echo "$version" | sed "s/.*SLES[\ -]\([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      elif test "$suse_open" ; then
        osver=`echo "$version" | sed "s/.*SUSE[\ -]\([^\.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      else
        osver=`echo "$version" | sed "s/.*[Ll][Ii][Nn][Uu][Xx]\ \(.*\)\ .*/\\1/"`
      fi
    elif test -r /etc/novell-release ; then
      ostype="suse"
      version=`grep -i linux /etc/novell-release`
      suse_ed=`echo "$version" | grep -i "Enterprise Desktop"`
      suse_es=`echo "$version" | grep -i "Enterprise Server"`
      suse_sles=`echo "$version" | grep -i "SLES"`
      if test "$suse_ed" ; then
        osver=`echo "$version" | sed "s/.*Desktop\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      elif test "$suse_es" ; then
        osver=`echo "$version" | sed "s/.*Server\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      elif test "$suse_sles" ; then
        osver=`echo "$version" | sed "s/.*SLES\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      else
        osver=`echo "$version" | sed "s/.*[Ll][Ii][Nn][Uu][Xx]\ \(.*\)\ .*/\\1/"`
      fi
    elif test -r /etc/sles-release ; then
      ostype="suse"
      version=`grep -i linux /etc/sles-release`
      suse_ed=`echo "$version" | grep -i "Enterprise Desktop"`
      suse_es=`echo "$version" | grep -i "Enterprise Server"`
      suse_sles=`echo "$version" | grep -i "SLES"`
      if test "$suse_ed" ; then
        osver=`echo "$version" | sed "s/.*Desktop\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      elif test "$suse_es" ; then
        osver=`echo "$version" | sed "s/.*Server\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      elif test "$suse_sles" ; then
        osver=`echo "$version" | sed "s/.*SLES\ \([^\ \.]*\)[\ \.].*/\\1/"`
        osver="SLE$osver"
      else
        osver=`echo "$version" | sed "s/.*[Ll][Ii][Nn][Uu][Xx]\ \(.*\)\ .*/\\1/"`
      fi
    elif test `cat /etc/issue | awk '{print $1}'` = "Ubuntu" ; then
      ostype=`cat /etc/issue | awk '{print $1}'`
      osver=`grep DISTRIB_RELEASE /etc/lsb-release | awk -F= '{print $2}'`
    elif test -r /etc/yellowdog-release ; then
      ostype="ydl"
      osver=`echo "$version" | sed "s/.*release\ \(.*\)\ .*/\\1/"`
    fi

    osmaj=`echo $osver | cut -d. -f1`
    osmin=`echo $osver | cut -d. -f2`
    if test "$osmaj" = "$osver" ; then
      osmaj=""
      osmin=""
    fi

    case `uname -m` in
      ppc64)  if test "$osver" = "SLES8" ; then
                arch="ppc32"
              else
                arch="ppc64"
              fi
              ;;
      x86_64) arch="x64"     ;;
      *86)    arch="x86"     ;;
      ppc)    arch="ppc32"   ;;
      ia64)   arch="ia64"    ;;
      s390)   arch="s390"    ;;
      *)      arch="unknown" ;;
    esac
    ;;

  SunOS)
    osver=`uname -r`
    osmaj=`echo $osver | cut -d. -f1`
    osmin=`echo $osver | cut -d. -f2`

    if test "$osmaj" -ge 5 ; then
      ostype="solaris"
      osmaj=`expr $osmaj - 3`
      osver="$osmaj.$osmin"

      if test -x /opt/SUNWspro/bin/fpversion ; then
        arch="sparc32"
      else
        arch="x64"
      fi
    fi
    ;;

  AIX)
    ostype="aix"
    osmaj=`uname -v`
    osmin=`uname -r`
    osver="$osmaj.$osmin"
    arch="ppc64"
    ;;

  *)
    osrep="unknown"
    ;;
esac

if test -z "$osrep" ; then
  if test -n "$arch" ; then
    osrep="$ostype-$osver-$arch"
  else
    osrep="$ostype-$osver"
  fi
fi


echo $osrep
