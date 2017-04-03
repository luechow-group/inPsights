#!/usr/bin/ruby

fmax = ARGV[0].to_i

if fmax < 3 then
  $stderr.puts "No terms for < 3rd degree"
  exit 1
end

total = 0

(3..fmax).each do |degree|
  terms = []
  puts "Degree #{degree}:"

  # free terms without r_ij
  (2..degree/2).each do |k|
    l = degree - k
    terms << "ri^#{k} rj^#{l} + ri^#{l} rj^#{k}"
  end

  # free terms with linear r_ij
  (2..(degree-1)/2).each do |k|
    l = degree - k
    m = degree - 1
    terms << "rij (ri^#{m} + rj^#{m} - ri^#{k} rj^#{l-1} - ri^#{l-1} rj^#{k})"
  end

  # free terms with higher order r_ij
  (2..degree-2).each do |k|
    l = degree - k
    terms << "rij^#{l} (ri^#{k} + rj^#{k})"
    (2..l/2).each do |m|
      terms << "rij^#{k} (ri^#{m} rj^#{l-m} + ri^#{l-m} rj^#{m})"
    end
  end

  # non-free terms
  nf = []
  l = degree - 1

  # non-free term without r_ij
  nf.push "ri rj^#{l} + ri^#{l} rj"

  # non-free term with linear r_ij
  nf.push "#{l == 2 ? "-0.5 " : "-"}rij (ri^#{l} + rj^#{l} - ri rj^#{l-1} - ri^#{l-1} rj)"

  # non-free term with higher order r_ij
  nf.push "rij^#{l} (ri + rj)"

  # non-free terms with products of r_i and r_j
  (2..degree-2).each do |m|
    k = l - m
    nf.push "#{k == 1 ? "0.5 " : ""}rij^#{m} (ri rj^#{k} + ri^#{k} rj)"
  end

  # generate free terms from non-free
  (0...nf.size-1).each do |t|
    terms << nf[t] + " - " + nf[-1]
  end

  terms.each do |term|
    puts "\t#{term}"
  end

  puts "\tTerms: #{terms.size}"
  total += terms.size
end

puts "Total: #{total}"
