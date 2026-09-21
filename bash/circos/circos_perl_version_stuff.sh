#!/bin/bash -l
module purge
module load Circos/0.69-10-GCCcore-13.3.0
perl -v | head -2        # confirm v5.38.2

# install cpanm into ~/perl5 using the correct Perl
curl -L https://cpanmin.us | perl - App::cpanminus --local-lib=$HOME/perl5

# put ~/perl5 on PERL5LIB / PATH for this shell
eval "$(perl -I $HOME/perl5/lib/perl5 -Mlocal::lib)"

# install the missing modules
cpanm Clone Config::General Font::TTF::Font Math::Bezier Params::Validate Math::Round List::MoreUtils Readonly Math::VecStat Regexp::Common Text::Format Set::IntSpan SVG Statistics::Basic

circos -modules

# git pull;circos -conf circos.conf