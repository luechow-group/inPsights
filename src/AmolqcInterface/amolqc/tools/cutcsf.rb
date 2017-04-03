#!/usr/bin/env ruby

if ARGV.size < 2 then
	puts <<EOD


    Usage: cutcsf.rb wfFile action cutoff EandSymlist


    "action" can be :

       "cielist" which print ci coefs and sum of orbital energies
            that form determinant for each csf.
       "csfelist" is same but it also prints determinants.
       "symlist"  prints out  symmetry in adition to csdelist
       "cicut" cut the expansion based on  abs(ci).
       "ecut" cut the expansion based on  some of energies.
       "convert" converts CSFS to determinants.
       "punchgms" punchs dat file from wf file in GAMESS format
       "punchci" punchs ci coeefs
       "settozero" set all ci coeeficients to zero
       "spcsf"    combine spin CSFs to build spin and spatial CSFs(only for diatomic molecules).

EOD
	exit
end
class Array
  def unordered_hash
    unless @_compare_o && @_compare_o == hash
      p = Hash.new(0)
      each{ |v| p[v] += 1 }
      @_compare_p = p.hash
      @_compare_o = hash
    end
    @_compare_p
  end
  def compare(b)
    unordered_hash == b.unordered_hash
  end
end

class DET
  attr_accessor :scpl, :detlist, :energy,:sym,:ci
end

class CSF
  attr_accessor :ci, :NumDets, :det
end
class ORB
  attr_accessor :energy, :sym
end
calcE=true
norm=0
orbs=[]
numcsf=0
wflines = IO.readlines ARGV[0]
action = ARGV[1]
if (action=="cielist" or action=="csfelist" or action=="symlist") then
  if ARGV.size < 3 then
      puts <<EOD

    for actions "cielist", "csfelist" and "symlist"
      you need 3 arguments.
      example:
        cutcsf.rb cielist wavefunction.wf csflist Elistfile.txt
EOD
      exit
  end
      eFile = IO.readlines(ARGV[2])
      eFile.each{|line|
        orb=ORB.new
        orb.energy= line.split.at(0).to_f
        orb.sym= line.split.at(1)
        orbs <<orb
        }
end
if (action=="cicut") then
  calcE=false
   if ARGV.size < 3 then
      puts <<EOD
      for action "cicut"
        you need 3 arguments.
        example:
        cutcsf.rb wavefunction.wf cicut 0.01
EOD
      exit
    end
    trs = ARGV[2].to_f
end
if (action=="ecut") then
   if ARGV.size < 4 then
      puts <<EOD
      for action "ecut"
        you need 4 arguments.
        example:
        cutcsf.rb wavefunction.wf ecut -3.58 Elistfile.txt
EOD
      exit
    end
    trs = ARGV[2].to_f
      eFile = IO.readlines(ARGV[3])
      eFile.each{|line|
        orb=ORB.new
        orb.energy= line.split.at(0).to_f
        orb.sym= line.split.at(1)
        orbs <<orb
        }
end
if (action=="convert" or action=="punchgms" or  action=="spcsf" or action=="settozero" or action=="punchci") then
  calcE=false
   if ARGV.size < 2 then
      puts <<EOD
      for actions "convert", "spcsf" and "settozero"
        you need 2 arguments.
        example:
        cutcsf.rb wavefunction.wf convert
EOD
      exit
    end
end

csfstart = wflines.find_index{|l| l.include? "$csfs" }
mostart  = wflines.find_index{|l| l.include? "$mos" }
ncsfs = wflines[csfstart + 1].to_i
crcsfLine = csfstart + 3
csfs = []

for i in 1..ncsfs
  ci, num_dets = wflines[crcsfLine-1].split

  csf = CSF.new
  csf.ci = ci.to_f
  csf.NumDets = num_dets = num_dets.to_i

  csf.det = wflines[crcsfLine ... crcsfLine + num_dets].map do |line|
    det_data = line.split

    det = DET.new
    det.scpl = det_data.shift.to_f
    det.detlist = det_data.map &:to_i
    det.energy = det.detlist.map{|i|orbs[i - 1].energy}.inject(:+) if calcE==true

    det
  end

  crcsfLine += num_dets + 1

  csfs << csf
end

if action == "symlist"
      csfno=1
     csfs.each{|csf|
      print "csfno: ", csfno,"\n"
      csfno+=1
        print "  ci coef = ", csf.ci,"   Sum of orbital energies  = ", csf.det[0].energy, "\n","\n"
        csf.det.each{|det| print "    ", det.scpl, "   "
         det.detlist.each{|list| print list,"(",orbs[list.to_i-1].sym, ") " }
         #det.detlist.each{|list| print orbs[list.to_i-1].sym," " }
            print   "\n","\n"}
        print "---\n"

     }

  end

if  action == "csfelist"
     csfs.each{|csf|
       puts [
         csf.ci,
         nil,
         csf.det.map{|det|
           "    #{det.scpl}   #{det.detlist.join " "}     #{det.energy}"
         },
         "---"
       ]
     }
end
if action == "ecut" then
  wflines.each_with_index{|line,index| puts line if index< csfstart }
	puts "$csfs"
	     csfs.each{|csf| if csf.det[0].energy < trs then
                         numcsf+=1
                         norm+=(csf.ci)**2
                        end
       }
       norm=Math.sqrt(norm)
       csfs.map{|csf| csf.ci=csf.ci/norm}
	     puts numcsf
       csfs.each{|csf|
       	if csf.det[0].energy < trs then
         puts [
           "#{csf.ci} " " #{csf.NumDets}",
           csf.det.map{|det|
             "  #{det.scpl}   #{det.detlist.join " "}"
           },
         ]
       end
        }
        puts "$end"
     $stderr.puts "CI coeefs renormalized after truncation! "
     $stderr.puts " initial norm = #{norm}"
end

if action == "cicut" then
    wflines.each_with_index{|line,index| puts line if index< csfstart }
		puts "$csfs"
		   csfs.each{|csf| if csf.ci.abs > trs then
                          numcsf+=1
                          norm+=(csf.ci)**2
                        end
          }
          norm=Math.sqrt(norm)
          csfs.map{|csf| csf.ci=csf.ci/norm }
		puts numcsf
       csfs.each{|csf|
       	if csf.ci.abs > trs/norm then
         puts [
           "#{csf.ci} " " #{csf.NumDets}",
           csf.det.map{|det|
             "  #{det.scpl}   #{det.detlist.join " "}"
           },
         ]
       end
        }
     puts "$end"
     $stderr.puts "CI coeefs renormalized after truncation! "
     $stderr.puts " initial norm = #{norm}"
end
if action == "punchci" then
    #wflines.each_with_index{|line,index| puts line if index< csfstart }
    #puts "ci coeffs"
       csfs.each{|csf|
                          norm+=(csf.ci)**2

          }
          norm=Math.sqrt(norm)
          csfs.map{|csf| csf.ci=csf.ci/norm }
       csfs.each{|csf|
         puts [
           "#{csf.ci} "
         ]
        }
     $stderr.puts "CI coeefs renormalized after truncation! "
     $stderr.puts " initial norm = #{norm}"
end

if action == "cielist"
     csfs.each{|csf|

       print  csf.ci,"   ", csf.det[0].energy, "\n"

     }
end

if action == "convert"
dets=[]
detscompact=[]
tmp=DET.new
csfs.each{|csf| csf.det.each{|det|
             tmp = det
             tmp.ci = det.scpl * csf.ci
                   dets << tmp
          }                  }

detscompact[0]=dets[0]
dets
for i in 1..dets.size-1
   considered =false
    for j in 0..detscompact.size-1
          if dets[i].detlist==detscompact[j].detlist then
            detscompact[j].ci = detscompact[j].ci + dets[i].ci
            considered=true
          end
    end
        if considered == false
            detscompact << dets[i]
        end
  end
  detscompact=detscompact.sort_by{|det| -(det.ci.abs)}
wflines.each_with_index{|line,index| puts line if index< csfstart }
puts "$dets"
puts detscompact.size
detscompact.each{|det|print det.ci," "
det.detlist.each{|list| print " ",list}
print "\n"

  }
  print "$end \n"

end

if action == "punchgms"

puts " $VEC"
wflines.each_with_index{|line,index| puts line if (index< csfstart-1 and index >mostart+2) }
puts " $END"

end



if action == "spcsf"
csfscomb=[]
wflines.each_with_index{|line,index| puts line if index< csfstart }
puts "$csfs"
       csfs.each{|csf|
                         norm+=(csf.ci)**2
       }
       norm=Math.sqrt(norm)
       csfs.map{|csf| csf.ci=csf.ci/norm}
     $stderr.puts "CI coeefs renormalized! "
     $stderr.puts " initial norm = #{norm}"

sw =false
cnt=-1
for i in 0..ncsfs-1
    if sw == true then
      sw = false
      next
    end
    if i<=cnt then
      next
    end
    rep=0
    for j in i+1..ncsfs-1
       if (csfs[i].ci.abs- csfs[j].ci.abs).abs<0.000001 and csfs[i].NumDets == csfs[j].NumDets then
         rep+=1
         idx=j
       end
    end

    if rep == 0 then
       csfscomb << csfs[i]
      elsif rep > 1
        $stderr.puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        $stderr.puts "!! To many CSFs with similar ci coeefs check CSF #{i+1} "
        $stderr.puts "!! The ci coeef is  #{csfs[i].ci} "
        $stderr.puts "--------------------------------------------------------"
        cnt=i+rep
        for l in i..i+rep
           csfscomb << csfs[l]
        end
      else
        csf=CSF.new
          csf.ci=csfs[i].ci
          csf.NumDets= csfs[i].NumDets*2
          d=[]
        csfs[i].det.each {|a| d << a }
        csfs[idx].det.each { |a| d << a }
        nrm=0
        d.each{|a|  nrm+= a.scpl**2}
        nrm=Math.sqrt(nrm)
        csf.ci=csf.ci*nrm
        d.map{|a| a.scpl=a.scpl/nrm}
        if csfs[i].ci * csfs[idx].ci < 0.0 then
            for l in csfs[i].NumDets..(csfs[i].NumDets*2)-1
              d[l].scpl=-d[l].scpl
            end
        end
         sw = true
         csf.det=d
         csfscomb << csf
       end

end
    puts csfscomb.size
       csfscomb.each{|csf|

         puts [
           "#{csf.ci} " " #{csf.NumDets}",
           csf.det.map{|det|
             "  #{det.scpl}   #{det.detlist.join " "}"
           },
         ]

        }
        puts "$end"

  $stderr.puts "!!CSFs should carfully investigated after combination!! "
  $stderr.puts "!!!!!!!!!!!no guarantee for correct resuts!!!!!!!!!!!!! "
end

if action == "settozero" then
    wflines.each_with_index{|line,index| puts line if index< csfstart }
    puts "$csfs"
    puts csfs.size
       csfs.each{|csf|

         puts [
           "0.000000 " " #{csf.NumDets}",
           csf.det.map{|det|
             "  #{det.scpl}   #{det.detlist.join " "}"
           },
         ]
        }
     puts "$end"
end
