# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021

require './periodic_table'

class Atom
  include PeriodicTable
  attr_reader :mass, :symbol, :color
  attr_accessor :id, :partition
  attr_writer :edges

  def initialize(id, symbol)
    @id = id
    @symbol = symbol
    @mass = PeriodicTable::ELEMENTS.index(@symbol) + 1
    @color = PeriodicTable::ELEMENT_COLORS[@symbol]
    @edges = []
    @partition = 0
  end

  def edges
    @edges.sort!.reverse! # return in descending order
  end
end

module Inchi

  def read_molfile(filename)
    molfile = File.read(filename) # reads entire file and closes it
    molfile.split("\n")
  end

  def create_molecule_array(molfile_lines)
    molecule = []
    atom_count, edge_count = molfile_lines[3].scan(/\d+/).map(&:to_i) # on 4th  line, 1st number is number of atoms, 2nd number is number of bonds.
    (4..atom_count + 3).each_with_index do |atom_index, i|
      molecule.push(Atom.new(i, molfile_lines[atom_index].split(' ')[3]))
    end
    (0..edge_count - 1).each do |edge_index|
      vertex1, vertex2 = parse_edge(molfile_lines[edge_index + 4 + atom_count])
      molecule[vertex1].edges.push(vertex2)    # add to the first atom of a bond
      molecule[vertex2].edges.push(vertex1)    # and to the second atom of the bond
    end
    molecule
  end

  def canonicalize_molecule(molecule, filename)
    filename = File.basename(filename, '.mol')

    print_molecule(molecule,
      "\nInitial data structure of #{filename}:")

    sorted_molecule = sort_atoms_by_mass(molecule)
    update_atom_ids(sorted_molecule) # indices need to be updated in order to be able to index molecule with atom id
    print_molecule(sorted_molecule,
      "\n#{filename} with atoms sorted by atomic mass (increasing):")

    partitioned_molecule = partition_molecule_by_mass(sorted_molecule)
    print_molecule(partitioned_molecule,
      "\n#{filename} partitioned by mass:")

    sorted_molecule
  end

  def write_ninchi_string(molecule)
    sum_formula = write_sum_formula_string(molecule)
    serialized_molecule = serialize_molecule(molecule)
    "nInChI=1S/#{sum_formula}/c#{serialized_molecule}"
  end

  def write_dot_file(molecule, filename)
    filename = File.basename(filename, '.mol')
    dotfile = "graph #{filename}\n{\n  bgcolor=grey\n"
    molecule.each_with_index do |atom, i|
      dotfile += "  #{i} [label=\"#{atom.symbol} #{i}\" color=#{atom.color},style=filled,shape=circle,fontname=Calibri];\n"
    end
    graph = compute_graph(molecule)
    graph.each { |line| dotfile += "  #{line[0]} -- #{line[1]} [color=black,style=bold];\n" if line[0] != line[1] }
    dotfile += "}\n"
  end

  # Helper methods ###################################################################################
  private

  def serialize_molecule(molecule)
    graph = compute_graph(molecule)
    inchi_string = ''
    graph.each do |line|
      inchi_string += "(#{line[0]}-#{line[1]})" if line[0] != line[1]
    end
    inchi_string
  end

  def compute_graph(molecule)
    graph = []
    molecule.each do |atom|
      atom.edges.each do |edge|
        graph.push([atom.id, edge].sort)
      end
    end
    graph.uniq.sort
  end

  def write_sum_formula_string(molecule)
    # Write sum formula in the order C > H > all other elements in alphabetic order.
    element_counts = compute_element_counts(molecule)
    element_counts.transform_values! { |v| v > 1 ? v : '' } # remove 1s since counts of 1 are implicit in sum formula
    sum_formula_string = ''
    sum_formula_string += "C#{element_counts['C']}" if element_counts.key?('C')
    sum_formula_string += "H#{element_counts['H']}" if element_counts.key?('H')
    element_counts.sort.to_h.each do |element, count|
      sum_formula_string += "#{element}#{count}" unless %w[C H].include?(element)
    end
    sum_formula_string
  end

  def compute_element_counts(molecule)
    # Compute hash table mapping element symbols to stoichiometric counts.
    unique_elements = molecule.map { |atom| atom.symbol }.uniq
    initial_counts = Array.new(unique_elements.length, 0)
    element_counts = unique_elements.zip(initial_counts).to_h
    molecule.each { |atom| element_counts[atom.symbol] += 1 }
    element_counts
  end

  def mass_edges(molecule, atom)
    # A Molecule class would be handy.
    atom.edges.map { |edge| molecule[edge].mass }.sort.reverse
  end

  def partition_molecule_by_mass(molecule)
    # Mutates molecule.
    current_partition = molecule[0].partition
    (0..molecule.size - 2).each do |i|
      j = i + 1
      fingerprint_atom_i = [molecule[i].mass] + mass_edges(molecule, molecule[i])
      fingerprint_atom_j = [molecule[j].mass] + mass_edges(molecule, molecule[j])
      current_partition += 1 if fingerprint_atom_i != fingerprint_atom_j
      molecule[j].partition = current_partition
    end
    molecule
  end

  def sort_atoms_by_mass(molecule)
    # Mutates molecule.
    molecule.sort do |atom_a, atom_b|
      fingerprint_atom_a = [atom_a.mass] + mass_edges(molecule, atom_a)
      fingerprint_atom_b = [atom_b.mass] + mass_edges(molecule, atom_b)
      fingerprint_atom_a <=> fingerprint_atom_b
    end
  end

  def update_atom_ids(molecule, random_indices=false, update_edges=true)
    # Mutates molecule.
    index_updates = compute_index_updates(molecule, random_indices)
    molecule.each do |atom|
      atom.id = index_updates[atom.id]
      atom.edges.map! { |edge| index_updates[edge] } if update_edges
    end
    molecule
  end

  def compute_index_updates(molecule, random_indices)
    current_indices = molecule.map(&:id)
    updated_indices = (0..molecule.length - 1).to_a
    updated_indices.shuffle! if random_indices
    current_indices.zip(updated_indices).to_h
  end

  def parse_edge(molfile_line)
    vertex1, vertex2, * = molfile_line.split(' ').map { |i| i.to_i - 1 }
    vertex1, vertex2 = vertex2, vertex1 if vertex1 > vertex2    # make sure first atom always has lower (not: higher?) index
    [vertex1, vertex2]
  end

  def print_molecule(molecule, caption)

    pad_width_id = molecule.map { |atom| atom.id.to_s.length }.max
    pad_width_mass = molecule.map { |atom| atom.mass.to_s.length }.max

    partitions = molecule.map { |atom| atom.partition.to_s }
    ids_masses = molecule.map { |atom| "(#{atom.id.to_s.ljust(pad_width_id)},#{atom.mass.to_s.rjust(pad_width_mass)})" }
    edge_ids_masses = molecule.map { |atom| atom.edges.zip(mass_edges(molecule, atom)).map { |id, mass| "(#{id.to_s.ljust(pad_width_id)},#{mass.to_s.rjust(pad_width_mass)})"} }
    edge_ids_masses.map! { |a| a.size == 1 ? a[0] : a.reduce('') { |s, si| s + "#{si}; " }.delete_suffix('; ') }

    col_headers = ['partition', '(index, mass)', '(index, mass) of connected atoms'.ljust(edge_ids_masses.map(&:length).max)]
    puts caption
    puts
    puts col_headers.reduce('|') { |header, col| header + " #{col} |" }
    puts
    (0..molecule.size - 1).each do |i|
      line = '|'
      line += " #{partitions[i].ljust(col_headers[0].length)} |"
      line += " #{ids_masses[i].ljust(col_headers[1].length)} |"
      line += " #{edge_ids_masses[i].ljust(col_headers[2].length)} |"
      puts line
    end
  end
end
