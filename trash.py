def make_umap_plot_from_data(sample, reference, umap_data_frame, close_up):
    """Create and save umap plot from UMAP data.

    Args:
        sample: sample data
        reference: reference data
        umap_data_frame: pandas data frame containing umap info. First
            row corresponds to sample.
        close_up: bool to indicate if only top matches should be plotted.
    """
    umap_sample = umap_data_frame.iloc[0]
    umap_title = f"UMAP for {sample.name} against {reference.name}, "\
        + f"{len(reference.annotated_specimens)} reference cases, "\
        + f"{len(sample.cpg_overlap)} CpGs"
    if close_up:
        umap_title = "Close-up " + umap_title
    umap_plot = px.scatter(
        x=umap_data_frame['x'],
        y=umap_data_frame['y'],
        labels={"x":"UMAP 0", "y":"UMAP 1", "color":"WHO class"},
        title=umap_title,
        color=umap_data_frame['methylation_class'],
        hover_name=umap_data_frame['id'],
        render_mode=config.PLOTLY_RENDER_MODE,
    )
    umap_plot.add_annotation(
        x=umap_sample['x'],
        y=umap_sample['y'],
        text=sample.name,
        showarrow=True,
        arrowhead=1,
    )
    umap_plot.update_yaxes(scaleanchor = "x", scaleratio = 1)
    # If close-up add hyperlinks for all references and draw circle
    if close_up:
        umap_plot.update_traces(marker=dict(size=5))
        # Add hyperlinks
        for _, row in umap_data_frame.iloc[1:].iterrows():
            umap_plot.add_annotation(
                x=row['x'],
                y=row['y'],
                text="<a href='" + config.CNV_URL_PREFIX + row['id']
                    + config.CNV_URL_SUFFIX
                    + "' target='_blank'>&nbsp;</a>",
                showarrow=False,
                arrowhead=1,
            )
        # Draw circle
        radius = umap_data_frame['distance'].iloc[-1]
        umap_plot.add_shape(
            type="circle",
            x0=umap_sample['x'] - radius,
            y0=umap_sample['y'] - radius,
            x1=umap_sample['x'] + radius,
            y1=umap_sample['y'] + radius,
            line_color="black",
            line_width=0.5,
        )
    path = umap_output_path(sample, reference, close_up)
    umap_plot.write_html(path['html'], config=dict({'scrollZoom': True}))
    umap_plot.write_image(path['png'])

def make_umap_report(sample, reference, close_up_data_frame):
    """Create PDF-compatible HTML table with proper cell padding, no borders
    etc."""
    # Generate Table
    html_table = f"""
        <table border=0>
        <tr>
            <td><b>methClass</b></td>
            <td><b>distance</b></td>
            <td><b>txt</b></td>
            <td><b>sentrix_ID</b></td>
        </tr>
        """
    for _, i in close_up_data_frame.iterrows():
        html_table += f"""
        <tr>
            <td>{i['methylation_class']}</td>
            <td>{i['distance']}</td>
            <td>{i['id']}</td>
            <td>{i['methylation_class']}</td>
        </tr>
        """

    # Generate PDF report
    html_report = f"""
        <html>
            <head>
                <title>{sample.name}</title>
            <body>
                <h1>{sample.name}</h1>{sample.name}
                <br>
                {html_table}
            </body>
        """
    file_path = os.path.join(
        config.NANODIP_REPORTS,
        f"{sample.name}_{reference.name}_NanoDiP_ranking.pdf",
    )
    convert_html_to_pdf(html_report, file_path)

    file_path = os.path.join(config.NANODIP_REPORTS, f"{sample.name}_cpgcount.txt")
    with open(file_path, "w") as f:
        f.write("%s" % len(sample.cpg_overlap))

def make_umap_plot(sample, reference):
    """Create html, png, csv files for UMAP methylation analysis.

    Args:
        sample: sample to analyse
        reference: reference data
    """
    import umap #TODO move to beginning

    logger.info(
        f"UMAP Plot initiated for {sample.name} and reference {reference.name}."
    )
    logger.info(f"Reference Annotation:\n{reference.annotation}")
    logger.info(f"Reference CpG Sites No:\n{len(reference.cpg_sites)}")
    logger.info(f"Reference Specimens No:\n{len(reference.specimens)}")
    logger.info(
        f"Reference Annotated specimens: {len(reference.annotated_specimens)}"
    )
    logger.info(f"Sample CpG Sites No:\n{len(sample.cpg_sites)}")
    logger.info(f"Sample CpG overlap No before:\n{sample.cpg_overlap}")

    # Calculate overlap of sample CpG's with reference CpG's (some probes have
    # been skipped from the reference set, e.g. sex chromosomes).
    sample.set_cpg_overlap(reference)
    logger.info(f"Sample read. CpG overlap No after:\n{len(sample.cpg_overlap)}")

    if sample.cpg_overlap:
        # Extract reference and sample methylation according to CpG overlap.
        reference_methylation = data.get_reference_methylation(sample,
                                                               reference)
        logger.info(f"""Reference methylation extracted:
                    {reference_methylation}""")
        sample_methylation = data.get_sample_methylation(sample, reference)
        logger.info(f"""Sample methylation extracted:
                    {sample_methylation}""")
        logger.info("UMAP algorithm initiated.")

        # Calculate UMAP Nx2 Matrix. Time intensive (~1min).
        umap_2d = umap.UMAP(verbose=True).fit_transform(
            np.vstack([sample_methylation, reference_methylation])
        )

        # Free memory
        del reference_methylation
        del sample_methylation

        logger.info("UMAP algorithm done.")

        umap_sample = umap_2d[0]
        umap_df = pd.DataFrame({
            'distance': [np.linalg.norm(z - umap_sample) for z in umap_2d],
            'methylation_class':  [sample.name] + reference.methylation_class,
            'id': [sample.name] + reference.specimen_ids,
            'x': umap_2d[:,0],
            'y': umap_2d[:,1],
        })

        file_path = os.path.join(config.NANODIP_REPORTS, f"{sample.name}_UMAP.csv")
        umap_df[['id', 'x', 'y']].to_csv(file_path, index=False)
        make_umap_plot_from_data(sample, reference, umap_df, close_up=False)
        logger.info("UMAP plot generated.")

        close_up_df = umap_df.sort_values(
            by='distance'
        )[:config.UMAP_PLOT_TOP_MATCHES + 1]
        make_umap_plot_from_data(sample, reference, close_up_df, close_up=True)
        logger.info("UMAP close-up plot generated.")

        make_umap_report(sample, reference, close_up_df)

    else:
        path = umap_output_path(sample, reference, close_up=True)
        with open(path['html'], 'w') as f:
            f.write("<html><body>No data to plot.</body></html>")

    logger.info("UMAP done.")

    @cherrypy.expose
    def ___umap(self, samp, ref):
        try:
            sample = data.SampleData(samp)
            reference = data.ReferenceData(ref)
        except FileNotFoundError:
            raise cherrypy.HTTPError(405, "URL not allowed")
        file_path = f"/reports/{samp}_{ref}_UMAP_all.html"
        print("FILE",file_path)
        #raise cherrypy.HTTPRedirect(file_path)
        if not os.path.exists(file_path):
            print("PATH does not exist***************")
            make_umap_plot(sample, reference)
        else:
            print("PATH exist******************")
            make_umap_plot(sample, reference)
        raise cherrypy.HTTPRedirect(file_path)

        t0=time.time()#TODO
        self.annotated_specimens_index = [self.specimens.index(a)
            for a in self.annotated_specimens]
        print("7", time.time()-t0);t0=time.time()
        self.annotated_specimens_index.sort()

        t0=time.time()#TODO
        self.methylation_class = []
        self.specimen_ids = []
        for i in self.annotated_specimens_index:
            s = self.specimens[i]
            self.specimen_ids.append(s)
            self.methylation_class.append(
                self.annotation.methylation_class[
                    self.annotation.id == s].values[0]
            )
        print("11", time.time()-t0);t0=time.time()

import jinja2

def url_for(static, filename):
    return f"{static}/{filename}"

url_for('static', filename='style.css')

def render_template(template_name, **context):
    loader = jinja2.FileSystemLoader("templates")
    template = jinja2.Environment(
        loader=loader).get_template(template_name)
    return template.render(context, url_for=url_for)

test_name = "FAQ17395_pass_barcode10_89406014_8"
test_dir = os.path.join(os.getcwd(), test_name)
sample_methylation = os.path.join(test_dir, test_name + "-freq.tsv")
output_overlap = os.path.join(test_dir, test_dir + "-methoverlap.tsv")
output_overlap_cnt = os.path.join(test_dir, test_dir + "-methoverlapcount.txt")


binary_search([1,3,4,8,11],-1,upper=True)
step_pos, copy_numbers, interval = [0, 1, 5, 6, 10], [0, 1, 1, 0], [5, 6]
median(step_pos, copy_numbers, interval)
read_intervals = read_positions
sp, cp = get_copy_numbers(genome.length, read_intervals)

l=100
reads = [[10,50],[40,60],[80,90]] # ([10, 40, 50, 60, 80, 90, 100], [0, 1, 2, 1, 0, 1, 0])
reads = [[0,70],[10,50],[40,60],[80,90]] #([10, 40, 50, 60, 70, 80, 90, 100], [1, 2, 3, 2, 1, 0, 1, 0])
reads = [[0,70],[10,50],[40,60],[80,90],[60,100]] #([10, 40, 50, 60, 70, 80, 90, 100], [1, 2, 3, 2, 2, 1, 2, 1])
reads = [[10,50],[50,60]] #([10, 50, 60, 100], [0, 1, 1, 0])
get_copy_numbers(l, reads)

def median(step_pos, copy_numbers, interval):
    """Calculate median from coverage step function."""
    if not len(interval) == 2:
        raise ValueError("invalid range")
    if not len(step_pos) == len(copy_numbers) + 1:
        raise ValueError("invalid lengths of step_pos/copy_numbers")
    left = binary_search(step_pos, interval[0])
    right = binary_search(step_pos, interval[1], upper=True)
    if left >= right:
        return 0
    delta_pos = np.diff(([max(interval[0], step_pos[left])]
                            + step_pos[left + 1:right]
                            + [min(interval[1], step_pos[right])]
                            ))
    multiplicity = list(zip(delta_pos, copy_numbers[left:right]))
    multiplicity.sort(key=lambda x: x[1])
    multiplicity_cum = np.cumsum([x[0] for x in multiplicity])
    total_len = multiplicity_cum[-1]
    i_med = binary_search(multiplicity_cum, total_len / 2, upper=True)
    if total_len % 2 == 1:
        return multiplicity[i_med][1]
    if total_len % 2 == 0:
        i_med_ = binary_search(multiplicity_cum, total_len / 2 + 1, upper=True)
        return (multiplicity[i_med][1] + multiplicity[i_med_][1]) / 2

# data
z=[[1, 20, 30],
  [20, 1, 60],
  [30, 60, 1]]

# every colorscale available to go.Heatmap() by default
scales =     ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
             'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
             'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
             'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
             'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
             'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
             'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
             'orrd', 'oryel', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg',
             'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn', 'puor',
             'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy',
             'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar', 'spectral',
             'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose',
             'tempo', 'temps', 'thermal', 'tropic', 'turbid', 'twilight',
             'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd']

# scales = scales[10:13]
colors =[]
buttons = []

# buttons, temporary figures and colorscales
for i, scale in enumerate(scales):
    colors.append(go.Figure(data=go.Heatmap(z=z,colorscale = scale)).data[0].colorscale)
    buttons.append(dict(method='restyle',
                        label=scale,
                        visible=True,
                        args=[{'colorscale':[colors[i]],
                                }, ],
                        )
                  )

# Initial figure:
fig = go.Figure(data=go.Heatmap(z=z, colorscale = scales[0]))

# some adjustments to the updatemenus
updatemenu = []
your_menu = dict()
updatemenu.append(your_menu)
updatemenu[0]['buttons'] = buttons
fig.update_layout(showlegend=False, updatemenus=updatemenu)

# f = fig.
fig.show()



figure(figsize=(20, 6), dpi=150)

plt.scatter(
    bin_midpoints,
    cnv,
    s=1,
    linewidths=5,
    c=copy_numbers,
    cmap=plt.cm.coolwarm_r,
    #vmin=cleanYLower,
    #vmax=cleanYUpper,
)
plt.show()

    def insert_(self, value):
        # Linear time complexity
        if not self.list or self.list[0] < value:
            self.list.insert(0, value)
        else:
            index = (len(self.list)
                    - next(i for i, x in enumerate(reversed(self.list))
                             if self.key(x) > self.key(value)))
            self.list.insert(index, value)

def make_binary_reference_data_(input_dir=REFERENCE_BETA_VALUES,
                               output_dir=REFERENCE_METHYLATION_DATA,
                               methylation_cutoff=METHYLATION_CUTOFF):
    """Create binary methylation files from raw reference data.

    Args:
        input_dir: Directory of reference data as float arrays-
            files.
        output_dir: Output directory containing binary array-files.
        methylation_cutoff: Empirical cutoff value for methylated
            (round to 1) and unmethylated (round to 0) CpGs.
    """

    specimens = [f for f in os.listdir(input_dir)
                 if f.endswith(".bin")]

    rounded_data = None
    for specimen in specimens:
        raw_data_path = os.path.join(input_dir, specimen)

        with open(raw_data_path, 'rb') as f:
            raw_data = np.fromfile(f, dtype=float)
            rounded_data = np.digitize(
                               raw_data,
                               bins=[methylation_cutoff]
                           ).astype(bool)


        output_path = os.path.join(output_dir, specimen)
        with open(output_path, 'wb') as f:
            f.write(rounded_data)

# Working from Home:
REFERENCE_BETA_VALUES = os.path.join(os.getcwd(), "data/reference/betaEPIC450Kmix_bin")
REFERENCE_METHYLATION_DATA = os.path.join(os.getcwd(), "data/reference/EPIC450K")
REFERENCE_ANNOTATIONS = os.path.join(os.getcwd(), "data/reference/reference_annotations")
NANODIP_OUTPUT = os.path.join(os.getcwd(), "data/nanodip_output")

annotated_reference_cases = []
not_annotated_reference_cases = []
for a in annotation["id"]:
    filename = a + "_betas_filtered.bin"
    path = os.path.join(REFERENCE_BETA_VALUES, filename)
    if os.path.exists(path):
        annotated_reference_cases.append(filename)
    else:
        not_annotated_reference_cases.append(filename)

analysis_reference = np.empty([len(annotated_reference_cases),
                               len(cpg_overlap)],
                              dtype=bool)


for i, case in enumerate(annotated_reference_cases):
    path = os.path.join(REFERENCE_BETA_VALUES, case)
    with open(path, "rb") as f:
        analysis_reference[i] = np.fromfile(f, dtype=bool)[cpg_overlap_index]

reference_methylation = np.fromfile(REFERENCE_METHYLATION, dtype=bool).reshape(shape)
analysis_reference = reference_methylation[grid]
del reference_methylation
