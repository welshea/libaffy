rm -rf .sconf_temp
rm config.log
rm ./libutils/config.log
find . -iname "*.scon*" -print0 | xargs -0 rm -rf
find . -iname "*.o"     -print0 | xargs -0 rm
find . -iname "*.a"     -print0 | xargs -0 rm
find . -iname "*.exe"   -print0 | xargs -0 rm
rm -rf bin/*
rm affy-apps/affydump/affydump
rm affy-apps/calvindump/calvindump
#rm affy-apps/datExtractor/datExtractor
rm affy-apps/findmedian/findmedian
rm affy-apps/iron/iron
rm affy-apps/iron_generic/iron_generic
rm affy-apps/mas5/mas5
rm affy-apps/pairgen/pairgen
rm affy-apps/rma/rma
