# moments_to_amplitudes.py
import ROOT
import os
import numpy as np
from typing import List, Dict, Tuple, Optional

class Moments2Amplitudes:
    """
    A class for analyzing moments-to-amplitudes conversion from ROOT files.
    Handles multiple fit results, chi-squared analysis, and plotting of amplitude vs mass.
    """
    
    def __init__(self, 
                 fit_result_dir: str,
                 output_dir: str,
                 filename_pattern: str = "Results*.root",
                 mass_name: str = "Mpi0eta",
                 n_mass_bins: int = 30,
                 mass_bin_width: float = 0.04,
                 first_mass_bin_center: float = 0.82):
        """
        Initialize the analyzer with directory paths and mass binning parameters.
        
        Args:
            fit_result_dir: Directory containing fit result files
            output_dir: Directory for output plots
            filename_pattern: Pattern for fit result files, use [] for file numbers
                            (e.g., "Results*.root" or "result_tree_BWa2_freeSWaves_ZeroPWaves_[]_nominal.root")
            mass_name: Name of mass variable in the tree
            n_mass_bins: Number of mass bins
            mass_bin_width: Width of each mass bin
            first_mass_bin_center: Center of first mass bin
        """
        self.fit_result_dir = fit_result_dir
        self.output_dir = output_dir
        self.filename_pattern = filename_pattern
        self.mass_name = mass_name
        self.n_mass_bins = n_mass_bins
        self.mass_bin_width = mass_bin_width
        self.first_mass_bin_center = first_mass_bin_center
        
        # Ensure output directory exists
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)

        # Physics notation dictionaries
        self.refl_dict = {"a": "+", "b": "-"}
        self.l_dict = {"0": "S", "1": "P", "2": "D", "3": "F", "4": "G"}
        
        # Storage for loaded data
        self.root_files = []
        self.file_paths = []
        self.tree_results = []
        self.tree_mass_dependent_params = []
        self.df_results = []
        self.chi2_values = {}
        self.lowest_chi2_files = []
        self.gr_amp_vs_mass = {}
        
        # Configure ROOT
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(1)
    
    def load_fit_files(self, file_pattern: str = None, max_files: int = 100) -> None:
        """
        Load multiple ROOT files based on a pattern.
        
        Args:
            file_pattern: Pattern with [] placeholder for file numbers. 
                         If None, uses self.filename_pattern from constructor
            max_files: Maximum number of files to load
        """
        if file_pattern is not None:
            self.filename_pattern = file_pattern
            
        pattern = self.filename_pattern
        file_names = [pattern.replace("[]", str(i)) for i in range(max_files)]
        file_paths = [os.path.join(self.fit_result_dir, name) for name in file_names]
        
        self.root_files = []
        self.file_paths = []
        self.chi2_values = {}
        for file_path in file_paths:
            if os.path.exists(file_path):
                root_file = ROOT.TFile(file_path, "READ")
                if not root_file or not root_file.IsOpen():
                    continue
                print(f"Loaded: {file_path}")
                chi2_tree = root_file.Get("chi2")
                if not chi2_tree:
                    print(f"No chi2 tree found in {file_path}")
                    root_file.Close()
                    continue
                chi2_tree.GetEntry(0)
                chi2_leaf = chi2_tree.GetLeaf("chi2")
                chi2_value = chi2_leaf.GetValue() if chi2_leaf else -1.0
                is_valid_leaf = chi2_tree.GetLeaf("isValid")
                is_valid_minimum = is_valid_leaf.GetValue() if is_valid_leaf else 1.0

                good_fit = (chi2_value > 0) and (is_valid_minimum == 1)
                if good_fit:
                    self.file_paths.append(file_path)
                    self.chi2_values[file_path] = chi2_value
                else:
                    print(f"Excluding file with invalid fit: {file_path} (chi2={chi2_value}, isValid={is_valid_minimum})")
                root_file.Close()

    # Backwards compatibility alias
    load_multiple_fit_files = load_fit_files

    def extract_chi2_values(self) -> Dict[str, float]:
        """
        Extract chi-squared values from all loaded ROOT files.
        
        Returns:
            Dictionary mapping file names to chi2 values
        """
        if not self.chi2_values:
            self.load_fit_files()
        
        return self.chi2_values
    
    def _get_chi2_from_file(self, root_file: ROOT.TFile) -> Optional[float]:
        """Extract chi2 value from a ROOT file."""
        chi2_tree = root_file.Get("chi2")
        if chi2_tree:
            chi2_tree.GetEntry(0)
            leaf = chi2_tree.GetLeaf("chi2")
            if leaf:
                return leaf.GetValue()
        return None
    
    def print_chi2_values(self, sort_by_value: bool = True, show_all: bool = True, n_display: int = None) -> None:
        """
        Print chi2 values in a formatted table.
        
        Args:
            sort_by_value: If True, sort by chi2 value (lowest first), else by filename
            show_all: If True, show all files, else show only n_display files
            n_display: Number of files to display (ignored if show_all=True)
        """
        if not self.chi2_values:
            self.extract_chi2_values()
        
        if not self.chi2_values:
            print("No chi2 values found!")
            return
        
        # Sort the values
        if sort_by_value:
            sorted_items = sorted(self.chi2_values.items(), key=lambda x: x[1])
        else:
            sorted_items = sorted(self.chi2_values.items(), key=lambda x: x[0])
        
        # Limit display if requested
        if not show_all and n_display is not None:
            sorted_items = sorted_items[:n_display]
        
        # Print header
        print("=" * 80)
        print(f"{'File Name':<50} {'Chi2 Value':>15} {'Rank':>10}")
        print("=" * 80)
        
        # Print values
        for rank, (file_name, chi2_value) in enumerate(sorted_items, 1):
            # Extract just the filename without path
            basename = os.path.basename(file_name)
            print(f"{basename:<50} {chi2_value:>15.4f} {rank:>10}")
        
        print("=" * 80)
        print(f"Total files: {len(self.chi2_values)}")
        if not show_all and n_display is not None:
            print(f"Showing: {min(n_display, len(sorted_items))} files")
        
        # Print statistics
        chi2_values_list = list(self.chi2_values.values())
        print(f"Min chi2: {min(chi2_values_list):.4f}")
        print(f"Max chi2: {max(chi2_values_list):.4f}")
        print(f"Mean chi2: {np.mean(chi2_values_list):.4f}")
        print(f"Std chi2: {np.std(chi2_values_list):.4f}")
        print("=" * 80)

    def select_best_fits(self, n_best: int = 20) -> List[Tuple[str, float]]:
        """
        Select the best fits based on lowest chi-squared values.
        
        Args:
            n_best: Number of best fits to select
            
        Returns:
            List of tuples (filename, chi2_value) sorted by chi2
        """
        if not self.chi2_values:
            self.extract_chi2_values()
        
        sorted_chi2 = sorted(self.chi2_values.items(), key=lambda x: x[1])
        self.lowest_chi2_files = sorted_chi2[:n_best]
        
        print(f"Selected {n_best} files with lowest chi2 values:")
        for file_name, chi2_value in self.lowest_chi2_files:
            print(f"{file_name}: {chi2_value:.2f}")
        
        return self.lowest_chi2_files

    def load_trees_from_best_fits(self, n_best: int = 20) -> None:
        """Load result and mass_dependent_params trees from the best fit files."""
        self.tree_results = []
        self.tree_mass_dependent_params = []
        self.root_files = []
        if (len(self.lowest_chi2_files) == 0): 
            self.select_best_fits(n_best)
        
        for file_name, _ in self.lowest_chi2_files:
            if os.path.exists(file_name):
                print(f"Loading tree from {file_name}")
                root_file = ROOT.TFile(file_name, "READ")
                if root_file.IsOpen():
                    tree = root_file.Get("result")
                    if tree:
                        entries_before = tree.GetEntries()
                        print(f"Original tree entries: {entries_before}")
                        # Clone the tree to avoid issues when file is closed
                        tree_clone = tree.CloneTree()
                        tree_clone.SetDirectory(0)  # Detach from file
                        entries_after = tree_clone.GetEntries()
                        print(f"Cloned tree entries: {entries_after}")
                        self.tree_results.append(tree_clone)
                    else:
                        print(f"Warning: No 'result' tree found in {file_name}")

                    tree_md = root_file.Get("mass_dependent_params")
                    if tree_md:
                        entries_before = tree_md.GetEntries()
                        print(f"Original mass_dependent_params tree entries: {entries_before}")
                        tree_clone_md = tree_md.CloneTree()
                        tree_clone_md.SetDirectory(0)
                        entries_after = tree_clone_md.GetEntries()
                        print(f"Cloned mass_dependent_params tree entries: {entries_after}")
                        self.tree_mass_dependent_params.append(tree_clone_md)

                    root_file.Close()
            else:
                print(f"Error: Could not open {file_name}")

    def cleanup_root_files(self) -> None:
        """Close all open ROOT files."""
        if hasattr(self, 'root_files'):
            for root_file in self.root_files:
                if root_file and root_file.IsOpen():
                    root_file.Close()
            self.root_files = []

    def load_df_result_from_best_fits(self, n_best_fits) -> None:
        """Load result dataframes from the best fit files."""
        self.df_results = []
        print("Loading dataframes from best fit files...")
        if (len(self.lowest_chi2_files) == 0): self.select_best_fits(n_best_fits)
        for file_name, _ in self.lowest_chi2_files:
            # print(f"Loading dataframe from {file_name}")
            df = ROOT.RDataFrame("result", file_name)
            self.df_results.append(df)

    def get_column_names(self, df_index: int = 0) -> List[str]:
        """Get column names from a ROOT dataframe."""
        if not self.df_results:
            self.load_df_result_from_best_fits(1)

        df = self.df_results[df_index]
        return list(df.GetColumnNames())
    
    def filter_amplitude_columns(self, col_names: List[str]) -> List[str]:
        """Filter column names to get only amplitude columns."""
        return [col for col in col_names if not col.startswith("aphi_") 
                and not col.startswith("bphi_") 
                and not col.endswith("_error")
                and not col.startswith("mass_bin") 
                and not col.startswith("H_")]
    
    def filter_moment_columns(self, col_names: List[str], moment_type: str = "H_") -> List[str]:
        """Filter column names to get moment columns."""
        return [col for col in col_names if col.startswith(moment_type)]

    def get_moment_name_list_L(self, L_values: List[int]) -> List[str]:
        """Filter column names to get moment columns based on specified alpha, L, M values."""
        moment_cols = self.get_column_names()
        filtered_moments = []
        for col in moment_cols:
            if col.startswith("H_"):
                parts = col.split("_")
                if len(parts) == 4:
                    alpha = parts[1]
                    L = int(parts[2])
                    M = parts[3]
                    if L in L_values:
                        filtered_moments.append(col)
        return filtered_moments

    def get_amplitude_column_names(self, df_index: int = 0) -> List[str]:
        """Get amplitude column names from a ROOT dataframe."""
        col_names = self.get_column_names(df_index)
        return self.filter_amplitude_columns(col_names)

    def get_moment_column_names(self, df_index: int = 0) -> List[str]:
        """Get moment column names from a ROOT dataframe."""
        col_names = self.get_column_names(df_index)
        return self.filter_moment_columns(col_names)

    def get_common_moment_names(self) -> Dict[str, str]:
        """
        Get dictionary of common moment names and their descriptions.
        
        Returns:
            Dictionary mapping moment names to descriptions
        """
        return {
            "H_0_0_0": "Moment (α=0, L=0, M=0)",
            "H_0_1_0": "Moment (α=0, L=1, M=0)",
            "H_0_1_1": "Moment (α=0, L=1, M=1)", 
            "H_0_2_0": "Moment (α=0, L=2, M=0)",
            "H_0_2_1": "Moment (α=0, L=2, M=1)",
            "H_0_2_2": "Moment (α=0, L=2, M=2)",
            "H_0_4_0": "Moment (α=0, L=4, M=0)",
            "H_0_4_1": "Moment (α=0, L=4, M=1)",
            "H_0_4_2": "Moment (α=0, L=4, M=2)",
            "H_0_4_3": "Moment (α=0, L=4, M=3)",
            "H_0_4_4": "Moment (α=0, L=4, M=4)",
            "H_1_1_0": "Moment (α=1, L=1, M=0)",
            "H_1_1_1": "Moment (α=1, L=1, M=1)",
            "H_1_2_0": "Moment (α=1, L=2, M=0)",
            "H_1_2_1": "Moment (α=1, L=2, M=1)",
            "H_1_2_2": "Moment (α=1, L=2, M=2)",
            "H_1_4_0": "Moment (α=1, L=4, M=0)",
            "H_1_4_1": "Moment (α=1, L=4, M=1)",
            "H_1_4_2": "Moment (α=1, L=4, M=2)",
            "H_1_4_3": "Moment (α=1, L=4, M=3)",
            "H_1_4_4": "Moment (α=1, L=4, M=4)"
        }
    
    def get_moments_L2(self, alpha_values: List[int] = [0, 1]) -> List[str]:
        """
        Get all L=2 moment names for specified alpha values.
        
        Args:
            alpha_values: List of alpha values to include
            
        Returns:
            List of moment names with L=2
        """
        moments_L2 = []
        for alpha in alpha_values:
            for M in range(3):  # M = 0, 1, 2 for L=2
                moments_L2.append(f"H_{alpha}_2_{M}")
        return moments_L2
    
    def get_moments_L4(self, alpha_values: List[int] = [0, 1]) -> List[str]:
        """
        Get all L=4 moment names for specified alpha values.
        
        Args:
            alpha_values: List of alpha values to include
            
        Returns:
            List of moment names with L=4
        """
        moments_L4 = []
        for alpha in alpha_values:
            for M in range(5):  # M = 0, 1, 2, 3, 4 for L=4
                moments_L4.append(f"H_{alpha}_4_{M}")
        return moments_L4
    
    def list_available_moments(self, moments_file_path: str) -> List[str]:
        """
        List all available moment graphs in a ROOT file.
        
        Args:
            moments_file_path: Path to ROOT file containing moment graphs
            
        Returns:
            List of moment graph names found in the file
        """
        if not os.path.exists(moments_file_path):
            print(f"Error: File not found: {moments_file_path}")
            return []
        
        root_file = ROOT.TFile(moments_file_path, "READ")
        if not root_file.IsOpen():
            print(f"Error: Could not open file: {moments_file_path}")
            return []
        
        moment_names = []
        try:
            for key in root_file.GetListOfKeys():
                hist_name = key.GetName()
                # Check for both graph patterns: H_ (TH1) and gr_unnormalized_moments_ (TGraphErrors)
                if hist_name.startswith("H_"):
                    graph = root_file.Get(hist_name)
                    if graph and isinstance(graph, ROOT.TH1):
                        moment_names.append(hist_name)
                elif hist_name.startswith("gr_unnormalized_moments_"):
                    graph = root_file.Get(hist_name)
                    if graph and isinstance(graph, ROOT.TGraphErrors):
                        moment_names.append(hist_name)
        finally:
            root_file.Close()
        
        return sorted(moment_names)

    def get_amplitude_sqr_vs_mass_graph(self, tree: ROOT.TTree, amp_name: str, unique_name: str) -> ROOT.TGraph:
        """
        Create a TGraph of amplitude squared vs mass.
        
        Args:
            tree: ROOT tree containing the data
            amp_name: Name of the amplitude column
            unique_name: Unique identifier for the graph name
            
        Returns:
            TGraph with amplitude squared vs mass
        """
        gr = ROOT.TGraphErrors()
        gr.SetName(f"{amp_name}_sqr_vs_mass_{unique_name}")
        gr.SetTitle(f"{amp_name} squared vs Mass ({unique_name})")
        
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("mass_bin", 1)
        # print(f"Setting branch status for amplitude: {amp_name}")
        tree.SetBranchStatus(ROOT.TString(amp_name).Data(), 1)
        # tree.SetBranchStatus(ROOT.TString(f"{amp_name}_error").Data(), 1)
        
        for entry in tree:
            mass = tree.GetLeaf("mass_bin").GetValue()
            amp_value = tree.GetLeaf(ROOT.TString(amp_name).Data()).GetValue()
            # amp_error = tree.GetLeaf(ROOT.TString(f"{amp_name}_error").Data()).GetValue()
            gr.SetPoint(gr.GetN(), mass, amp_value * amp_value)
            # gr.SetPointError(gr.GetN() - 1, 0, 2 * amp_value * amp_error)
        
        return gr
    
    def get_amplitude_sqr_vs_mass_hist(self, tree: ROOT.TTree, amp_name: str, unique_name: str) -> ROOT.TH1D:
        """
        Create a histogram of amplitude squared vs mass.
        
        Args:
            tree: ROOT tree containing the data
            amp_name: Name of the amplitude column
            unique_name: Unique identifier for the histogram name
        Returns:
            TH1D histogram with amplitude squared vs mass
        """
        h1 = ROOT.TH1D(
            f"h1_{amp_name}_sqr_vs_mass_{unique_name}", 
            f"{amp_name} squared vs Mass ({unique_name})", 
            self.n_mass_bins,
            self.first_mass_bin_center - self.mass_bin_width/2,
            self.first_mass_bin_center + (self.n_mass_bins - 1) * self.mass_bin_width + self.mass_bin_width/2
        )
        
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("mass_bin", 1)
        tree.SetBranchStatus(ROOT.TString(amp_name).Data(), 1)
        # tree.SetBranchStatus(ROOT.TString(f"{amp_name}_error").Data(), 1)
        for entry in tree:
            mass = tree.GetLeaf("mass_bin").GetValue()
            amp_value = tree.GetLeaf(ROOT.TString(amp_name).Data()).GetValue()
            # amp_error = tree.GetLeaf(ROOT.TString(f"{amp_name}_error").Data()).GetValue()
            h1.SetBinContent(h1.FindBin(mass), amp_value * amp_value)
            # h1.SetBinError(h1.FindBin(mass), 2 * amp_value * amp_error)
        
        return h1
    
    def get_moment_vs_mass_graph(self, tree: ROOT.TTree, moment_name: str) -> ROOT.TGraph:
        """Create a TGraph of moment vs mass."""
        gr = ROOT.TGraph()
        
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("mass_bin", 1)
        tree.SetBranchStatus(moment_name, 1)
        
        # Check if branches exist
        mass_branch = tree.GetBranch("mass_bin")
        moment_branch = tree.GetBranch(moment_name)
        
        if not mass_branch:
            print(f"Warning: mass_bin branch not found in tree")
            return gr
        if not moment_branch:
            print(f"Warning: {moment_name} branch not found in tree")
            return gr
        
        mass_bins = self._get_mass_bins()
        min_mass = self.first_mass_bin_center - self.mass_bin_width / 2.0
        max_mass = self.first_mass_bin_center + (self.n_mass_bins - 0.5) * self.mass_bin_width
        
        for entry in range(tree.GetEntries()):
            tree.GetEntry(entry)
            mass_bin_value = tree.GetLeaf("mass_bin").GetValue()
            moment_value = tree.GetLeaf(moment_name).GetValue()
            if min_mass <= mass_bin_value <= max_mass or mass_bin_value > 0:
                gr.SetPoint(gr.GetN(), mass_bin_value, moment_value/2)
        
        # Return a clone to avoid ROOT object management issues
        gr_clone = gr.Clone(f"{moment_name}_vs_mass_clone")
        return gr_clone
    
    def _get_mass_bins(self) -> List[float]:
        """Get mass bin centers."""
        return [self.first_mass_bin_center + self.mass_bin_width * i 
                for i in range(self.n_mass_bins)]
    
    def generate_amplitude_graphs_for_all_fits(self, unique_name: str = None) -> None:
        """Generate amplitude graphs for all loaded fits."""
        if not self.tree_results:
            self.load_trees_from_best_fits()
        
        col_names = self.get_column_names()
        amp_col_names = self.filter_amplitude_columns(col_names)
        
        self.gr_amp_vs_mass = {}
        for itree, tree in enumerate(self.tree_results):
            for amp_name in amp_col_names:
                print(f"Processing amplitude: {amp_name} for tree {itree}")
                gr = self.get_amplitude_sqr_vs_mass_graph(tree, amp_name, unique_name if unique_name else None)
                self.gr_amp_vs_mass[f"{amp_name}_{itree}"] = gr
    
    def plot_unnormalized_moments(self, moments_file_path: str, 
                                 moment_name: str, 
                                 output_filename: str,
                                 draw_best_fits: bool = True) -> None:
        """
        Plot unnormalized moments from graph ROOT files with optional best fit overlays.
        
        Args:
            moments_file_path: Path to ROOT file containing unnormalized moment graphs
            moment_name: Exact name of moment to plot. Format: H_α_L_M where:
                        - α: alpha parameter (0, 1, 2, ...)
                        - L: orbital angular momentum quantum number (0, 1, 2, 4, ...)
                        - M: magnetic quantum number (0, 1, 2, ..., L)
                        Examples: "H_0_0_0", "H_0_2_1", "H_1_4_3", "H_0_1_0"
                        Special cases: 
                        - "all" to plot all moments found in the file
                        - "L=2" to plot all moments with L=2 (any α, M)
                        - "L=4" to plot all moments with L=4 (any α, M)
            output_filename: Name of the output PDF file
            draw_best_fits: If True, overlay best fit results on the graphs
        """
        if not os.path.exists(moments_file_path):
            print(f"Error: File not found: {moments_file_path}")
            return
        
        root_file = ROOT.TFile(moments_file_path, "READ")
        if not root_file.IsOpen():
            print(f"Error: Could not open file: {moments_file_path}")
            return
        
        # Load trees from best fits if we want to draw them and haven't loaded them yet
        if draw_best_fits and not self.tree_results and self.lowest_chi2_files:
            self.load_trees_from_best_fits()
        
        try:
            # Check for special plotting modes
            if moment_name.lower() == "all":
                # Plot all moments - look for both H_ and gr_unnormalized_moments_ formats
                moment_graphs = []
                for key in root_file.GetListOfKeys():
                    hist_name = key.GetName()
                    print(f"Checking graph: {hist_name}")
                    if hist_name.startswith("H_"):
                        hist = root_file.Get(hist_name)
                        if hist and isinstance(hist, ROOT.TH1):
                            moment_graphs.append((hist_name, hist))
                    elif hist_name.startswith("gr_unnormalized_moments_"):
                        # Parse moment name to extract α, L, M values from gr_unnormalized_moments_α_L_M_vs_mass
                        parts = hist_name.split("_")
                        if len(parts) >= 6:  # gr_unnormalized_moments_α_L_M_vs_mass
                            try:
                                hist_alpha = int(parts[3])  # α is the 4th part (index 3)
                                hist_L = int(parts[4])      # L is the 5th part (index 4)  
                                hist_M = int(parts[5])      # M is the 6th part (index 5)
                                graph = root_file.Get(hist_name)
                                if graph and isinstance(graph, ROOT.TGraphErrors):
                                    moment_graphs.append((hist_name, graph))
                            except (ValueError, IndexError):
                                continue  # Skip malformed moment names
                
                canvas_title = "All Unnormalized Moments"
                
            elif moment_name.lower().startswith("l="):
                # Plot all moments with specific L value
                L_value = moment_name.lower().replace("l=", "")
                try:
                    L_int = int(L_value)
                except ValueError:
                    print(f"Error: Invalid L value '{L_value}'. Must be an integer.")
                    return
                
                moment_graphs = []
                for key in root_file.GetListOfKeys():
                    hist_name = key.GetName()
                    print(f"Checking graph: {hist_name}")
                    if hist_name.startswith("gr_unnormalized_moments_"):
                        # Parse moment name to extract L value from gr_unnormalized_moments_α_L_M_vs_mass
                        parts = hist_name.split("_")
                        if len(parts) >= 6:  # gr_unnormalized_moments_α_L_M_vs_mass
                            try:
                                hist_alpha = int(parts[3])  # α is the 4th part (index 3)
                                hist_L = int(parts[4])      # L is the 5th part (index 4)  
                                hist_M = int(parts[5])      # M is the 6th part (index 5)
                                if hist_L == L_int:
                                    graph = root_file.Get(hist_name)
                                    if graph and isinstance(graph, ROOT.TGraphErrors):
                                        moment_graphs.append((hist_name, graph))
                            except (ValueError, IndexError):
                                continue  # Skip malformed moment names
                
                canvas_title = f"All L={L_int} Moments"
                
            elif moment_name.startswith("H_"):
                # Plot single moment - look for exact match first
                alpha, L, M = moment_name.split("_")[1:4]
                moment_name = f"gr_unnormalized_moments_{alpha}_{L}_{M}_vs_mass"
                graph = root_file.Get(moment_name)
                if not graph:
                    print(f"Error: Graph '{moment_name}' not found in file")
                    print("Available moment graphs:")
                    available = self.list_available_moments(moments_file_path)
                    for avail_moment in available:
                        print(f"  - {avail_moment}")
                    return
                                   
                if not (isinstance(graph, ROOT.TH1) or isinstance(graph, ROOT.TGraphErrors)):
                    print(f"Error: Object '{moment_name}' is not a graph or histogram")
                    return
                
                canvas = ROOT.TCanvas("canvas_moment", f"Unnormalized Moment H_{alpha}_{L}_{M}", 
                                    800, 600)
                self._plot_single_moment_with_fits(graph, moment_name, draw_best_fits)
                
                canvas.Update()
                output_path = os.path.join(self.output_dir, output_filename)
                canvas.Print(output_path)
                print(f"Saved unnormalized moments plot to: {output_path}")
                return
            
            # Handle multi-moment plotting (all or L-specific)
            if not moment_graphs:
                if moment_name.lower() == "all":
                    print("No moment graphs found with 'H_' prefix")
                else:
                    print(f"No moments found with L={L_int}")
                return
            
            # Calculate canvas dimensions
            n_moments = len(moment_graphs)
            n_cols = min(4, n_moments)  # Max 4 columns
            n_rows = (n_moments + n_cols - 1) // n_cols  # Ceiling division
            
            canvas_width = n_cols * 400
            canvas_height = n_rows * 300
            canvas = ROOT.TCanvas("canvas_moments", canvas_title, 
                                canvas_width, canvas_height)
            canvas.Divide(n_cols, n_rows)
            
            for i, (hist_name, graph) in enumerate(moment_graphs):
                canvas.cd(i + 1)
                self._plot_single_moment_with_fits(graph, hist_name, draw_best_fits)
            
            canvas.Update()
            output_path = os.path.join(self.output_dir, output_filename)
            canvas.Print(output_path)
            print(f"Saved {canvas_title} plot to: {output_path}")
            
        finally:
            root_file.Close()
    
    def _plot_single_moment_with_fits(self, graph, moment_name: str, draw_best_fits: bool = True) -> None:
        """
        Helper method to plot a single moment graph with optional best fit overlays.
        
        Args:
            graph: ROOT graph object (TH1 or TGraphErrors)
            moment_name: Name of the moment for labeling
            draw_best_fits: If True, overlay best fit results
        """
        # Handle both TH1 and TGraphErrors
        if isinstance(graph, ROOT.TGraphErrors):
            # Clone graph to avoid issues with file closure
            graph_clone = graph.Clone(f"{graph.GetName()}_clone")
            
            # Set graph style
            graph_clone.SetLineColor(ROOT.kBlack)
            graph_clone.SetLineWidth(2)
            graph_clone.SetMarkerStyle(20)
            graph_clone.SetMarkerSize(0.8)
            graph_clone.SetMarkerColor(ROOT.kBlack)
            
            # Set axis labels
            graph_clone.GetXaxis().SetTitle("M_{#eta#pi^{0}} (GeV/c^{2})")
            graph_clone.GetYaxis().SetTitle(f"Unnormalized moment value")
            graph_clone.GetXaxis().SetTitleSize(0.05)
            graph_clone.GetYaxis().SetTitleSize(0.05)

            graph_clone.SetTitle(f"H_{{{moment_name.split('_')[3]}}}({moment_name.split('_')[4]},{moment_name.split('_')[5]})")
            
            # Draw graph
            graph_clone.Draw("APE")
            
        elif isinstance(graph, ROOT.TH1):
            # Clone histogram to avoid issues with file closure
            graph_clone = graph.Clone(f"{graph.GetName()}_clone")
            
            # Set histogram style
            graph_clone.SetLineColor(ROOT.kBlack)
            graph_clone.SetLineWidth(2)
            graph_clone.SetFillColor(ROOT.kCyan)
            graph_clone.SetFillStyle(3004)  # Hatched fill
            graph_clone.SetMarkerStyle(20)
            graph_clone.SetMarkerSize(0.8)
            graph_clone.SetMarkerColor(ROOT.kBlack)
            
            # Set axis labels
            graph_clone.GetXaxis().SetTitle("M_{#eta#pi^{0}} (GeV/c^{2})")
            graph_clone.GetYaxis().SetTitle(f"Unnormalized {moment_name}")
            graph_clone.GetXaxis().SetTitleSize(0.05)
            graph_clone.GetYaxis().SetTitleSize(0.05)
            
            # Draw histogram
            graph_clone.Draw("HIST E")
        else:
            print(f"Error: Unsupported object type for {moment_name}")
            return
        
        # Draw best fit results if requested and available
        if draw_best_fits and self.tree_results:
            col_names = self.get_column_names()
            moment_cols = self.filter_moment_columns(col_names, moment_type="H_")

            # convert moment_name to H_alpha_L_M format
            parts = moment_name.split("_")
            alpha = int(parts[3])
            L = int(parts[4])
            M = int(parts[5])
            moment_name = f"H_{alpha}_{L}_{M}"
            
            if moment_name in moment_cols:
                # Create graphs for each best fit result
                fit_graphs = []
                for i, tree in enumerate(self.tree_results):
                    gr = self.get_moment_vs_mass_graph(tree, moment_name)
                    if gr.GetN() > 0:  # Only add if graph has points
                        # Set different colors for different fits
                        gr.SetLineColorAlpha(ROOT.kRed, 0.3)
                        gr.SetLineWidth(2)
                        
                        # Draw the graph
                        gr.Draw("C same")
                        fit_graphs.append(gr)
                
                if fit_graphs:
                    print(f"Drew {len(fit_graphs)} best fit overlays for {moment_name}")
            else:
                print(f"Warning: Moment {moment_name} not found in fit results")
                print(f"Available moment columns: {moment_cols}")
        
        # Set margins
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetRightMargin(0.05)

    def _plot_single_moment(self, hist: ROOT.TH1, moment_name: str) -> None:
        """
        Helper method to plot a single moment histogram.
        
        Args:
            hist: ROOT histogram object
            moment_name: Name of the moment for labeling
        """
        # Clone histogram to avoid issues with file closure
        hist_clone = hist.Clone(f"{hist.GetName()}_clone")
        
        # Set histogram style
        hist_clone.SetLineColor(ROOT.kBlue)
        hist_clone.SetLineWidth(2)
        hist_clone.SetFillColor(ROOT.kBlue)
        hist_clone.SetFillStyle(3004)  # Hatched fill
        
        # Set axis labels
        hist_clone.GetXaxis().SetTitle("M_{#eta#pi^{0}} (GeV/c^{2})")
        hist_clone.GetYaxis().SetTitle(f"Unnormalized {moment_name}")
        hist_clone.GetXaxis().SetTitleSize(0.05)
        hist_clone.GetYaxis().SetTitleSize(0.05)
        
        # Draw histogram
        hist_clone.Draw("HIST")
        
        # Add text label
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.06)
        text.SetTextFont(42)
        text.DrawLatex(0.15, 0.85, moment_name)
        
        # Set margins
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetRightMargin(0.05)

    def get_single_unnormalized_moment_graph(self, moments_file_path: str, moment_name: str, unique_name: str = None) -> ROOT.TGraphErrors:
        """
        Get a single unnormalized moment graph from a ROOT file.
        
        Args:
            moments_file_path: Path to ROOT file containing moment graphs
            moment_name: Exact name of moment to retrieve
        Returns:
            TGraphErrors object of the requested moment
        """
        if not os.path.exists(moments_file_path):
            print(f"Error: File not found: {moments_file_path}")
            return None
        
        root_file = ROOT.TFile(moments_file_path, "READ")
        if not root_file.IsOpen():
            print(f"Error: Could not open file: {moments_file_path}")
            return None
        
        try:
            alpha, L, M = moment_name.split("_")[1:4]
            moment_name = "gr_unnormalized_moments_" + alpha + "_" + L + "_" + M + "_vs_mass"
            graph = root_file.Get(moment_name)
            if not graph or not isinstance(graph, ROOT.TGraphErrors):
                print(f"Error: Graph '{moment_name}' not found or is not a TGraphErrors in file")
                return None
            
            # Clone graph to avoid issues with file closure
            graph_clone = graph.Clone(f"{unique_name if unique_name else moment_name}_clone")
            return graph_clone
        finally:
            root_file.Close()
    
    def get_single_unnormalized_moment_hist(self, moments_file_path: str, moment_name: str, unique_name: str) -> ROOT.TH1D:
        """
        Get a single unnormalized moment histogram from a ROOT file.
        
        Args:
            moments_file_path: Path to ROOT file containing moment histograms
            moment_name: Exact name of moment to retrieve
            unique_name: Unique name for the output histogram
        Returns:
            TH1D object of the requested moment
        """
        graph = self.get_single_unnormalized_moment_graph(moments_file_path, moment_name)
        if not graph:
            return None
        
        # Convert TGraphErrors to TH1D
        n_points = graph.GetN()
        h1 = ROOT.TH1D(unique_name, f"{moment_name}", n_points, 
                       graph.GetX()[0] - self.mass_bin_width/2,
                       graph.GetX()[n_points-1] + self.mass_bin_width/2)
        
        for i in range(n_points):
            y = graph.GetPointY(i)
            h1.SetBinContent(i + 1, y)
            h1.SetBinError(i + 1, graph.GetErrorY(i))

        return h1

    def get_chi2_distribution_hist(self, n_bins: int, x_min: Optional[float] = None, x_max: Optional[float] = None, unique_name: Optional[str] = None) -> ROOT.TH1D:
        """
        Get the chi-squared distribution histogram.
        
        Returns:
            TH1D histogram of chi-squared values
        """
        if not self.chi2_values:
            self.extract_chi2_values()
        
        chi2_list = list(self.chi2_values.values())
        if x_min is None:
            x_min = min(chi2_list) * 0.95
        if x_max is None:
            x_max = max(chi2_list) * 1.05
        
        h1_chi2 = ROOT.TH1D(unique_name if unique_name!=None else "h1_chi2", "Chi2 Values Distribution", 
                           n_bins, x_min, x_max)
        
        for chi2_value in chi2_list:
            h1_chi2.Fill(chi2_value)
        
        return h1_chi2.Clone(unique_name if unique_name!=None else "h1_chi2_clone")

    def plot_chi2_distribution(self, output_filename: str, 
                              n_bins: int = 100, 
                              x_min: Optional[float] = None, 
                              x_max: Optional[float] = None) -> None:
        """
        Plot the distribution of chi-squared values.
        
        Args:
            output_filename: Name of the output PDF file
            n_bins: Number of histogram bins
            x_min, x_max: Histogram range (auto-determined if None)
        """
        if not self.chi2_values:
            self.extract_chi2_values()
        
        chi2_list = list(self.chi2_values.values())
        if x_min is None:
            x_min = min(chi2_list) * 0.95
        if x_max is None:
            x_max = max(chi2_list) * 1.05
        
        canvas = ROOT.TCanvas("canvas_chi2", "Chi2 Values", 1200, 800)
        h1_chi2 = ROOT.TH1D("h1_chi2", "Chi2 Values Distribution", 
                           n_bins, x_min, x_max)
        
        for chi2_value in chi2_list:
            h1_chi2.Fill(chi2_value)
        
        h1_chi2.SetFillColor(ROOT.kBlue)
        h1_chi2.SetLineColor(ROOT.kBlue)
        h1_chi2.GetXaxis().SetTitle("#chi^{2} value")
        h1_chi2.GetYaxis().SetTitle("Counts")
        h1_chi2.Draw("HIST")
        
        canvas.Update()
        canvas.Print(os.path.join(self.output_dir, output_filename))

    def get_H000_moment_graph(self, moments_file_path: str, unique_name: str) -> ROOT.TGraphErrors:
        """
        Get the H_0_0_0 moment graph from a ROOT file.
        
        Args:
            moments_file_path: Path to ROOT file containing moment graphs
        Returns:
            TGraphErrors object of the H_0_0_0 moment
        """
        return self.get_single_unnormalized_moment_graph(moments_file_path, "H_0_0_0", unique_name=f"gr_H_0_0_0_{unique_name}")
    
    def get_H000_moment_hist(self, moments_file_path: str, unique_name: str) -> ROOT.TH1D:
        """
        Get the H_0_0_0 moment histogram from a ROOT file.
        
        Args:
            moments_file_path: Path to ROOT file containing moment graphs
        Returns:
            TH1D object of the H_0_0_0 moment
        """
        return self.get_single_unnormalized_moment_hist(moments_file_path, "H_0_0_0", unique_name=f"h1_H_0_0_0_{unique_name}")
    
    def plot_amplitude_comparison(self, wave_type: str, reference_graph: ROOT.TGraph,
                                 output_filename: str, unique_name: str = None, use_histogram: bool = False, scale=1.0) -> None:
        """
        Plot amplitude comparison for a specific wave type.
        
        Args:
            wave_type: Wave type to plot (e.g., "a_2", "b_2", "ab_0")
            reference_graph: Reference graph (e.g., H_0_0_0 moment)
            output_filename: Output PDF filename
            unique_name: Unique name for the plot
            use_histogram: Whether to use histogram instead of graphs
        """
        col_names = self.get_column_names()
        amp_col_names = self.filter_amplitude_columns(col_names)

        if wave_type == "all":
            selected_names = [col for col in amp_col_names if not col.endswith("error")]
        elif wave_type == "ab_0":
            selected_names = [col for col in amp_col_names 
                            if (col.startswith("a_0") or col.startswith("b_0")) and not col.endswith("error")]
        else:
            selected_names = [col for col in amp_col_names if col.startswith(wave_type)]
        
        print(f"Selected {len(selected_names)} amplitudes for wave type: {wave_type}")
        
        if not selected_names:
            print(f"No amplitudes found for wave type: {wave_type}")
            return
        
        if (wave_type == "all"):
            canvas = ROOT.TCanvas(f"canvas_{wave_type}_{unique_name if unique_name else 'comparison'}", 
                            f"{wave_type} amplitudes vs Mass", 
                            800, 600)
            
            reference_graph.Draw("AP")
            x_min = self.first_mass_bin_center - self.mass_bin_width/2
            x_max = self.first_mass_bin_center + (self.n_mass_bins - 1) * self.mass_bin_width + self.mass_bin_width/2
            reference_graph.GetXaxis().SetRangeUser(x_min, x_max)
            reference_graph.SetMinimum(0.0)
            reference_graph.SetLineWidth(2)
            
            legend = ROOT.TLegend(0.75, 0.7, 1.0, 0.9)
            legend.AddEntry(reference_graph, "H_{0}(0,0)", "lpe")
            legend.SetTextSize(0.06)
            legend.SetBorderSize(0)
            legend.SetFillColorAlpha(0, 0.0)

            hist = self.get_amplitude_sqr_vs_mass_hist(self.tree_results[0], selected_names[0], unique_name if unique_name else None)

            for j, name in enumerate(selected_names):
                if j == 0:
                    continue  # Already drawn the first one
                
                hist.Add(self.get_amplitude_sqr_vs_mass_hist(self.tree_results[0], name, unique_name if unique_name else None))

            hist.SetLineColor(ROOT.kRed)
            hist.SetLineWidth(2)
            hist.Draw("HIST same")

            legend = ROOT.TLegend(0.75, 0.7, 1.0, 0.9)
            legend.AddEntry(reference_graph, "H_{0}(0,0)", "lpe")
            legend.SetTextSize(0.06)
            legend.SetBorderSize(0)
            legend.SetFillColorAlpha(0, 0.0)
            legend.AddEntry(hist, "#Sigma |l|^{2}", "l")
            legend.DrawClone()              

        else:
            n_plots = len(selected_names)
            if n_plots == 1:
                canvas_width = 800
            elif n_plots == 2:
                canvas_width = 1400
            elif n_plots == 3:
                canvas_width = 1800
            elif n_plots == 4:
                canvas_width = 2400
            elif n_plots == 5:
                canvas_width = 2800
            else:  # 6 or more
                canvas_width = n_plots * 500
            
            canvas = ROOT.TCanvas(f"canvas_{wave_type}_{unique_name if unique_name else ''}", 
                            f"{wave_type} amplitudes vs Mass", 
                            canvas_width, 600)
            canvas.Divide(len(selected_names), 1)
            for j, name in enumerate(selected_names):
                canvas.cd(j + 1)

                # Draw reference graph
                reference_graph.Draw("AP")
                x_min = self.first_mass_bin_center - self.mass_bin_width/2
                x_max = self.first_mass_bin_center + (self.n_mass_bins - 1) * self.mass_bin_width + self.mass_bin_width/2
                reference_graph.GetXaxis().SetRangeUser(x_min, x_max)
                reference_graph.SetMinimum(0.0)
                reference_graph.SetLineWidth(2)

                # Setup legend and pad margins
                legend = ROOT.TLegend(0.75, 0.7, 1.0, 0.9)
                legend.AddEntry(reference_graph, "H_{0}(0,0)", "lpe")
                legend.SetTextSize(0.06)
                legend.SetBorderSize(0)
                legend.SetFillColorAlpha(0, 0.0)

                ROOT.gPad.SetLeftMargin(0.15)
                ROOT.gPad.SetRightMargin(0.0)

                # Parse amplitude name for labeling
                print(f"Plotting amplitude: {name}")
                refl, l, m = "", "", ""
                if (len(name.split("_")) == 3):
                    parts = name.split("_")
                    refl = parts[0]  # a or b
                    l = parts[1]     # L value
                    m = parts[2]     # M value            

                # Plot amplitude data
                if use_histogram and len(self.tree_results) > 0:
                    hist = self.get_amplitude_sqr_vs_mass_hist(self.tree_results[0], name, unique_name=f"h1_{name}_{unique_name if unique_name else None}")
                    hist.Scale(scale)
                    hist.GetXaxis().SetTitle("M_{#eta#pi^{0}} (GeV/c^{2})")
                    hist.GetYaxis().SetTitle(f"|{self.l_dict[l]}^{{{self.refl_dict[refl]}}}_{{{m}}}|^{{{2}}}")
                    hist.SetLineColor(ROOT.kGray)
                    hist.SetLineWidth(2)
                    hist.SetFillColorAlpha(ROOT.kGray, 0.5)
                    hist.DrawClone("HIST same")
                    if (scale == 1.0):
                        legend.AddEntry(hist, f"{self.l_dict[l]}^{{{self.refl_dict[refl]}}}_{{{m}}}", "lf")
                    else:
                        legend.AddEntry(hist, f"{self.l_dict[l]}^{{{self.refl_dict[refl]}}}_{{{m}}} x {scale}", "lf")
                else:
                    for i in range(len(self.tree_results)):
                        if f"{name}_{i}" in self.gr_amp_vs_mass:
                            gr = self.gr_amp_vs_mass[f"{name}_{i}"]
                            gr.Scale(scale)
                            gr.SetLineColorAlpha(12, 0.2)
                            gr.SetLineWidth(2)
                            gr.GetXaxis().SetTitle("M_{#eta#pi^{0}} (GeV/c^{2})")
                            gr.GetYaxis().SetTitle(f"|{self.l_dict[l]}^{{{self.refl_dict[refl]}}}_{{{m}}}|^{{{2}}}")
                            gr.DrawClone("C same")
                            if i == 0:
                                if (scale == 1.0):
                                    legend.AddEntry(gr, f"{self.l_dict[l]}^{{{self.refl_dict[refl]}}}_{{{m}}}", "l")
                                else:
                                    legend.AddEntry(gr, f"{self.l_dict[l]}^{{{self.refl_dict[refl]}}}_{{{m}}} x {scale}", "l")

                legend.DrawClone()

        canvas.Draw()
        canvas.Print(os.path.join(self.output_dir, output_filename))

    def get_mass_dependent_parameters(self, 
                                     arg1 = None, 
                                     arg2 = None, 
                                     parameters: List[str] = None, 
                                     lowest_chi2_tree_index: int = 0) -> Dict[str, float]:
        """
        Get mass-dependent parameters from a specific tree index.
        Supports both signatures:
          get_mass_dependent_parameters(parameters, lowest_chi2_tree_index=0)
          get_mass_dependent_parameters(lowest_chi2_tree_index, parameters)
        
        Args:
            parameters: List of parameter names to extract
            lowest_chi2_tree_index: Index of the tree corresponding to the lowest chi-squared fit
        Returns:
            Dictionary of parameter names and their corresponding values
        """
        # Disambiguate arguments
        target_params = parameters
        tree_idx = lowest_chi2_tree_index

        if isinstance(arg1, (list, tuple)):
            target_params = list(arg1)
            if isinstance(arg2, int):
                tree_idx = arg2
        elif isinstance(arg1, int):
            tree_idx = arg1
            if isinstance(arg2, (list, tuple)):
                target_params = list(arg2)
        elif arg2 is not None and isinstance(arg2, (list, tuple)):
            target_params = list(arg2)

        if target_params is None:
            target_params = []

        if not self.tree_mass_dependent_params and not self.tree_results:
            self.load_trees_from_best_fits()
        
        tree_source = self.tree_mass_dependent_params if self.tree_mass_dependent_params else self.tree_results
        if not tree_source:
            print("Error: No trees available to extract mass-dependent parameters")
            return {}

        if tree_idx < 0 or tree_idx >= len(tree_source):
            print(f"Error: Invalid tree index {tree_idx}")
            return {}
        
        tree = tree_source[tree_idx]
        tree.SetBranchStatus("*", 0)
        for param in target_params:
            tree.SetBranchStatus(param, 1)
        
        param_values = {}
        tree.GetEntry(0)
        for param in target_params:
            leaf = tree.GetLeaf(param)
            if leaf:
                param_values[param] = leaf.GetValue()
            else:
                print(f"Warning: Parameter '{param}' not found in tree branches")
        
        return param_values

    # Compatibility alias
    GetMassDependentParameters = get_mass_dependent_parameters
