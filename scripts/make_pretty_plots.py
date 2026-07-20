"""
make_pretty_plots.py — Regenerate all project plots with improved styling.
Reads cleaned data directly; does NOT re-run notebooks.
"""
import os, math, warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.dates as mdates
import matplotlib as mpl
import seaborn as sns
from matplotlib.dates import YearLocator, DateFormatter
from scipy.stats import pearsonr
from shapely.geometry import Polygon
from shapely import wkt
import contextily as ctx

warnings.filterwarnings("ignore")

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLEAN = os.path.join(P, "data", "processed")
DATA  = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")

# ── Global rcParams ──────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.size':          11,
    'axes.titlesize':     13,
    'axes.labelsize':     11,
    'xtick.labelsize':    9,
    'ytick.labelsize':    9,
    'legend.fontsize':    9,
    'figure.dpi':         200,
    'savefig.dpi':        200,
    'savefig.bbox':       'tight',
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.grid':          False,
})

# ── Cohesive colour palette ──────────────────────────────────────────────
C_BLUE    = '#2166AC'   # gun homicide primary / "before"
C_RED     = '#D6604D'   # ShotSpotter primary
C_ORANGE  = '#F4A582'   # "during" period
C_BG_BEFORE = '#DEEBF7' # light blue wash
C_BG_DURING = '#FEE0D2' # light peach wash

SELECTED = ['Austin', 'Humboldt Park', 'North Lawndale',
            'Englewood', 'West Englewood', 'South Shore']

# ── Load data ────────────────────────────────────────────────────────────
print("Loading data …")
crime_six = pd.read_csv(os.path.join(CLEAN, "Homicide_six_hoods.csv"))
crime_six["Date"] = pd.to_datetime(crime_six["Date"])

ss_six = pd.read_csv(os.path.join(CLEAN, "Shotspotter_six_hoods.csv"))
ss_six["Date"] = pd.to_datetime(ss_six["Date"])

# Recompute shot_dates exactly as the notebook does
westside = ['Austin', 'Humboldt Park', 'North Lawndale']
southside = ['Englewood', 'West Englewood', 'South Shore']
ss_west = ss_six[ss_six['Community'].isin(westside)]
ss_south = ss_six[ss_six['Community'].isin(southside)]
neighborhoods_map = {h: ss_west for h in westside}
neighborhoods_map.update({h: ss_south for h in southside})

shot_dates = {}
for hood, df in neighborhoods_map.items():
    hd = df[df['Community'] == hood].sort_values('Date')
    first = hd.iloc[0]['Date']
    last  = hd.iloc[-1]['Date']
    dur   = last - first
    shot_dates[hood] = {
        'Start Date':      first - dur,
        'First Shot Date': first,
        'Last Shot Date':  last,
        'Duration': f"{dur.days // 365} years, "
                    f"{(dur.days % 365) // 30} months, "
                    f"{(dur.days % 365) % 30} days",
    }

def savefig(fig, name):
    p = os.path.join(PLOTS, name)
    fig.savefig(p, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    sz = os.path.getsize(p) / 1024
    print(f"  OK: {name}  ({sz:.0f} KB)")


# =====================================================================
# 1. homicide_comparison.png  — before / during bar chart
# =====================================================================
def plot_bar_chart():
    print("\n1. homicide_comparison.png")
    fig, axes = plt.subplots(2, 3, figsize=(14, 10))
    axes = axes.flatten()

    # compute counts
    counts = {}
    gmax = 0
    for hood in SELECTED:
        d = shot_dates[hood]
        hd = crime_six[crime_six['Community'] == hood]
        bef = hd[(hd['Date'] >= d['Start Date']) & (hd['Date'] < d['First Shot Date'])]
        dur = hd[(hd['Date'] >= d['First Shot Date']) & (hd['Date'] <= d['Last Shot Date'])]
        counts[hood] = (bef.shape[0], dur.shape[0])
        gmax = max(gmax, bef.shape[0], dur.shape[0])

    for idx, hood in enumerate(SELECTED):
        ax = axes[idx]
        bc, dc = counts[hood]
        bars = ax.bar(['Before', 'During'], [bc, dc],
                       color=[C_BLUE, C_ORANGE], width=0.45,
                       edgecolor='black', linewidth=0.5)
        ax.set_ylim(0, gmax * 1.12)
        ax.set_ylabel('Gun Homicides', fontsize=9)
        ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)
        ax.set_axisbelow(True)

        ax.set_title(hood, fontsize=12, fontweight='bold', pad=8)

        d = shot_dates[hood]
        txt = (f"Before: {d['Start Date'].strftime('%Y-%m-%d')} – "
               f"{d['First Shot Date'].strftime('%Y-%m-%d')}\n"
               f"During: {d['First Shot Date'].strftime('%Y-%m-%d')} – "
               f"{d['Last Shot Date'].strftime('%Y-%m-%d')}\n"
               f"Duration: {d['Duration']}")
        ax.text(0.5, -0.22, txt, fontsize=7, color='#555555',
                ha='center', va='center', transform=ax.transAxes)

        for bar, val in zip(bars, [bc, dc]):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + gmax * 0.015,
                    str(val), ha='center', va='bottom', fontsize=9, fontweight='bold')

    bpatch = mpatches.Patch(color=C_BLUE,   label='Before ShotSpotter')
    opatch = mpatches.Patch(color=C_ORANGE, label='During ShotSpotter')
    fig.legend(handles=[bpatch, opatch], loc='lower center', ncol=2,
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle('Gun Homicide Counts Before vs. During ShotSpotter',
                 fontsize=15, fontweight='bold', y=0.98)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig, "homicide_comparison.png")


# =====================================================================
# 2. Monthly_Gun_Homicide_ShotSpotter.png — 6-panel time series
# =====================================================================
def plot_monthly_timeseries():
    print("\n2. Monthly_Gun_Homicide_ShotSpotter.png")
    fig, axes = plt.subplots(6, 1, figsize=(18, 22))

    for idx, (hood, dates) in enumerate(shot_dates.items()):
        ax1 = axes[idx]
        sd, fsd, lsd = dates['Start Date'], dates['First Shot Date'], dates['Last Shot Date']

        # Homicides
        hd = crime_six[crime_six['Community'] == hood]
        mh = hd.groupby(hd['Date'].dt.to_period('M')).size().reset_index(name='cnt')
        mh['Date'] = mh['Date'].dt.to_timestamp()

        # ShotSpotter
        sp = ss_six[ss_six['Community'] == hood]
        ms = sp.groupby(sp['Date'].dt.to_period('M')).size().reset_index(name='cnt')
        ms['Date'] = ms['Date'].dt.to_timestamp()

        # Shading
        ax1.axvspan(mh['Date'].min(), fsd, alpha=0.15, color=C_BG_BEFORE, zorder=0)
        ax1.axvspan(fsd, mh['Date'].max(), alpha=0.15, color=C_BG_DURING, zorder=0)

        # Primary: homicides
        ax1.plot(mh['Date'], mh['cnt'], marker='o', markersize=2.5,
                 linewidth=1, color=C_BLUE, label='Monthly Gun Homicides', zorder=3)
        ax1.set_ylabel('Gun Homicides', color=C_BLUE, fontsize=9)
        ax1.tick_params(axis='y', labelcolor=C_BLUE, labelsize=8)

        # Secondary: ShotSpotter
        ax2 = ax1.twinx()
        ax2.spines['right'].set_visible(True)
        ax2.plot(ms['Date'], ms['cnt'], marker='x', markersize=2.5,
                 linewidth=0.8, alpha=0.7, color=C_RED,
                 label='Monthly ShotSpotter Alerts', zorder=2)
        ax2.set_ylabel('ShotSpotter Alerts', color=C_RED, fontsize=9)
        ax2.tick_params(axis='y', labelcolor=C_RED, labelsize=8)

        # Trend lines
        for subset, per in [(mh[mh['Date'] < fsd], 'Before'),
                            (mh[mh['Date'] >= fsd], 'During')]:
            if len(subset) >= 3:
                x = mdates.date2num(subset['Date']); y = subset['cnt']
                if np.any(y):
                    z = np.polyfit(x, y, 1); p = np.poly1d(z)
                    ax1.plot(subset['Date'], p(x), ':', color='#4A1486',
                             linewidth=1.2, alpha=0.6,
                             label=f'{per} Trend (Homicides)')
        for subset, per in [(ms[ms['Date'] < fsd], 'Before'),
                            (ms[ms['Date'] >= fsd], 'During')]:
            if len(subset) >= 3:
                x = mdates.date2num(subset['Date']); y = subset['cnt']
                if np.any(y):
                    z = np.polyfit(x, y, 1); p = np.poly1d(z)
                    ax2.plot(subset['Date'], p(x), ':', color='#B2182B',
                             linewidth=1.2, alpha=0.6,
                             label=f'{per} Trend (Shots)')

        ax1.axvline(fsd, color='#B2182B', ls='--', lw=1, zorder=4,
                    label='ShotSpotter Installation')

        # Axis limits
        if hood in ['Austin', 'North Lawndale']:
            ax1.set_ylim(0, 16)
        else:
            ax1.set_ylim(0, 10)
        ax2.set_ylim(0, 350)

        # X-axis: year ticks only, horizontal
        ax1.xaxis.set_major_locator(YearLocator(1))
        ax1.xaxis.set_major_formatter(DateFormatter('%Y'))
        ax1.tick_params(axis='x', labelsize=8, rotation=0)
        ax1.xaxis.set_minor_locator(plt.NullLocator())

        ax1.set_title(f'{hood}', fontsize=12, fontweight='bold', pad=6)
        ax1.text(0.5, -0.10,
                 f'ShotSpotter Installation: {fsd.strftime("%Y-%m-%d")}',
                 transform=ax1.transAxes, ha='center', fontsize=7, color='#666666')

        # Legend — upper left, small, semi-transparent bg
        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax1.legend(h1 + h2, l1 + l2, loc='upper left', fontsize=6,
                   frameon=True, framealpha=0.85, edgecolor='#cccccc',
                   ncol=2)

    fig.suptitle(
        'Monthly Gun Homicides & ShotSpotter Alerts\n'
        'with Before / During Trend Lines',
        fontsize=15, fontweight='bold', y=0.995)
    fig.subplots_adjust(hspace=0.35)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    savefig(fig, "monthly_timeseries.png")


# =====================================================================
# 3. Monthly mean counts — dual axis
# =====================================================================
def plot_monthly_mean():
    print("\n3. Monthly_Gun_Homicide_and_ShotSpotter_Mean_Counts_2017_2023.png")
    hom_post = crime_six[crime_six['Date'] >= '2017-01-13'].copy()
    mh = (hom_post.groupby(hom_post['Date'].dt.to_period('M'))
          .size().reset_index(name='cnt'))
    mh['Month'] = mh['Date'].dt.month
    hom_mean = mh.groupby('Month')['cnt'].mean().round().reset_index()

    ms = (ss_six.groupby(ss_six['Date'].dt.to_period('M'))
          .size().reset_index(name='cnt'))
    ms['Month'] = ms['Date'].dt.month
    ss_mean = ms.groupby('Month')['cnt'].mean().round().reset_index()

    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax1.plot(hom_mean['Month'], hom_mean['cnt'], 'o-',
             color=C_BLUE, lw=1.5, markersize=6, label='Gun Homicide Mean')
    ax1.set_xlabel('Month')
    ax1.set_ylabel('Mean Monthly Gun Homicides', color=C_BLUE)
    ax1.tick_params(axis='y', labelcolor=C_BLUE)
    ax1.set_xticks(range(1, 13))
    ax1.set_xticklabels(months, rotation=0)
    ax1.yaxis.grid(True, alpha=0.2, linewidth=0.5)
    ax1.set_axisbelow(True)
    ax1.set_ylim(hom_mean['cnt'].min() - 2, hom_mean['cnt'].max() + 2)

    ax2 = ax1.twinx()
    ax2.spines['right'].set_visible(True)
    ax2.plot(ss_mean['Month'], ss_mean['cnt'], 'o-',
             color=C_RED, lw=1.5, markersize=6, label='ShotSpotter Mean')
    ax2.set_ylabel('Mean Monthly ShotSpotter Alerts', color=C_RED)
    ax2.tick_params(axis='y', labelcolor=C_RED)
    ss_lo = ss_mean['cnt'].min() - 20
    ss_hi = ss_mean['cnt'].max() + 20
    ax2.set_ylim(ss_lo, ss_hi)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc='upper left', frameon=True,
               framealpha=0.9, edgecolor='#cccccc')

    ax1.set_title('Monthly Gun Homicide & ShotSpotter Mean Counts (2017–2023)',
                  fontsize=13, fontweight='bold', pad=10)
    fig.tight_layout()
    savefig(fig, "monthly_means.png")


# =====================================================================
# 4. Correlation scatter — 6 subplots
# =====================================================================
def plot_correlation():
    print("\n4. Correlation scatter plots")
    hom_post = crime_six[crime_six['Date'] >= '2017-01-13'].copy()
    hom_post['Year'] = hom_post['Date'].dt.year
    hom_post['Month'] = hom_post['Date'].dt.month

    md = hom_post.groupby(['Community', 'Year', 'Month']).size().reset_index(name='Homicide Count')
    sd = ss_six.groupby(['Community', 'Year', 'Month']).size().reset_index(name='ShotSpotter Alerts')
    corr = pd.merge(md, sd, on=['Community', 'Year', 'Month'], how='outer')
    corr['Homicide Count']     = corr['Homicide Count'].fillna(0)
    corr['ShotSpotter Alerts'] = corr['ShotSpotter Alerts'].fillna(0)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)
    axes = axes.flatten()

    for idx, hood in enumerate(SELECTED):
        ax = axes[idx]
        cd = corr[corr['Community'] == hood]
        r, pv = pearsonr(cd['Homicide Count'], cd['ShotSpotter Alerts'])

        sns.regplot(x='ShotSpotter Alerts', y='Homicide Count', data=cd, ax=ax,
                    scatter_kws={'s': 25, 'alpha': 0.4, 'color': C_BLUE,
                                 'edgecolors': 'none'},
                    line_kws={'color': C_RED, 'lw': 1.5},
                    ci=95,
                    color=C_RED)
        # Patch the confidence band to lighter red
        for coll in ax.collections:
            if hasattr(coll, 'get_facecolor'):
                fc = coll.get_facecolor()
                if len(fc) > 0 and fc[0][3] < 1:  # it's the band
                    coll.set_facecolor((*mpl.colors.to_rgb('#FDDBC7'), 0.25))

        # Proper minus sign
        r_str = f"{r:.2f}".replace("-", "\N{MINUS SIGN}")
        ax.set_title(f"{hood} (r = {r_str}, p = {pv:.2f})", fontsize=11, fontweight='bold')
        ax.set_xlim(0, 350)
        ax.set_ylim(0, 10)
        ax.yaxis.grid(True, alpha=0.15, linewidth=0.4)
        ax.xaxis.grid(True, alpha=0.15, linewidth=0.4)
        ax.set_axisbelow(True)
        ax.set_xlabel('Monthly ShotSpotter Alerts', fontsize=9)
        ax.set_ylabel('Gun Homicide Count', fontsize=9)

    fig.suptitle(
        'Correlation between ShotSpotter Alerts and\n'
        'Gun Homicide Counts by Community (2017–2023)',
        fontsize=14, fontweight='bold', y=1.01)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    savefig(fig, "correlation_by_community.png")


# =====================================================================
# 5. Hourly dual-axis
# =====================================================================
def plot_hourly():
    print("\n5. Hourly dual-axis plot")
    hom_al = crime_six[crime_six['Date'] >= '2017-01-13']
    hh = hom_al.groupby('Hour').size().reset_index(name='cnt')
    sh = ss_six.groupby('Hour').size().reset_index(name='cnt')

    fig, ax1 = plt.subplots(figsize=(12, 5.5))
    ax1.plot(sh['Hour'], sh['cnt'], 'o-', color=C_RED, markersize=5,
             lw=1.3, label='ShotSpotter Alerts')
    ax1.set_xticks(range(0, 24))
    ax1.set_xlabel('Hour of the Day')
    ax1.set_ylabel('ShotSpotter Alerts', color=C_RED)
    ax1.tick_params(axis='y', labelcolor=C_RED)
    ax1.yaxis.grid(True, alpha=0.2, linewidth=0.4)
    ax1.set_axisbelow(True)

    ax2 = ax1.twinx()
    ax2.spines['right'].set_visible(True)
    ax2.plot(hh['Hour'], hh['cnt'], 'o-', color=C_BLUE, markersize=5,
             lw=1.3, label='Gun Homicide Count')
    ax2.set_ylabel('Gun Homicide Count', color=C_BLUE)
    ax2.tick_params(axis='y', labelcolor=C_BLUE)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    # upper center: both series occupy the upper corners (alerts peak at hour 0-2,
    # homicides at hour 19-23), so the top-middle band is the only clear region.
    ax1.legend(h1 + h2, l1 + l2, loc='upper center', frameon=True,
               framealpha=0.9, edgecolor='#cccccc')

    ax1.set_title(
        'ShotSpotter Alerts & Gun Homicide Counts by Hour of the Day (2017–2023)',
        fontsize=13, fontweight='bold', pad=10)
    fig.tight_layout()
    savefig(fig, "hourly_pattern.png")


# =====================================================================
# 6. Murder-map choropleths (pre / post / full period)
# =====================================================================
def plot_murder_maps():
    print("\n6. Murder maps (pre-2017 / post-2017 / full period)")
    gdf = gpd.read_file(os.path.join(DATA, "CommAreas_20250408.geojson"))
    streets_raw = pd.read_csv(os.path.join(DATA, "transportation_20250415.csv"))
    homicide = pd.read_csv(os.path.join(CLEAN, "Homicide_six_hoods.csv"))

    # Pre-parse streets once
    streets_raw["geometry"] = streets_raw["the_geom"].apply(wkt.loads)
    streets_gdf = gpd.GeoDataFrame(streets_raw, geometry="geometry", crs="EPSG:4326")
    streets_gdf = streets_gdf.to_crs(gdf.crs)

    def make_map(hom_data, start, end, filename):
        filt = hom_data[(hom_data['Date'] >= start) & (hom_data['Date'] <= end)]
        ann = filt.groupby('Community').size().reset_index(name='Homicide Count')
        ann['Rank'] = ann['Homicide Count'].rank(method='first', ascending=False).astype(int)
        sel = dict(zip(ann['Community'].str.title(), ann['Rank']))

        hi = gdf[gdf['community'].str.title().isin(sel.keys())]
        xmin, ymin, xmax, ymax = gdf.total_bounds

        fig, ax = plt.subplots(figsize=(6, 8), dpi=200)
        ax.set_facecolor('#f5f5f5')
        ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)

        gdf.plot(ax=ax, color='#e8e8e8', edgecolor='#999999', linewidth=0.4,
                 alpha=0.8, zorder=1)
        hi.plot(ax=ax, color='#ff6347', edgecolor='#4d4d4d', alpha=0.7,
                zorder=2)
        streets_gdf.plot(ax=ax, color='#888888', linewidth=0.15, alpha=0.3,
                         zorder=3)

        for _, row in hi.iterrows():
            c = row.geometry.centroid
            num = sel.get(row['community'].title())
            ax.text(c.x, c.y, str(num), fontsize=13, fontweight='bold',
                    color='black', ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.65,
                              edgecolor='none'),
                    zorder=5)

        ax.set_title(f"Chicago Most Violent Neighborhoods\n{start} – {end}",
                     fontsize=12, fontweight='bold', pad=8)
        ax.axis('off')

        # North arrow
        ax.annotate('N', xy=(0.97, 0.95), xycoords='axes fraction',
                    fontsize=11, fontweight='bold', ha='center')
        ax.annotate('', xy=(0.97, 0.94), xycoords='axes fraction',
                    xytext=(0.97, 0.89), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.2))

        # Scale bar (correct metric)
        m_per_deg = 111320 * math.cos(math.radians(41.85))
        bar_deg = 500 / m_per_deg
        bx0 = xmin + (xmax - xmin) * 0.1
        by0 = ymin + (ymax - ymin) * 0.05
        ax.plot([bx0, bx0 + bar_deg], [by0, by0], 'k-', lw=2, zorder=6)
        ax.text(bx0 + bar_deg / 2, by0 + (ymax - ymin) * 0.012,
                '500 m', fontsize=8, ha='center', zorder=6)

        # Legend
        sn = sorted(sel.items(), key=lambda x: x[1])
        legend_text = "\n".join(f"{n}: {name}" for name, n in sn)
        ax.annotate(legend_text, xy=(0.01, 0.02), xycoords='figure fraction',
                    fontsize=8, ha='left', va='bottom',
                    fontfamily='monospace')

        savefig(fig, filename)

    make_map(homicide, "2009-01-03", "2017-01-13", "murder_map_pre_2017.png")
    make_map(homicide, "2017-01-13", "2023-12-31", "murder_map_post_2017.png")
    make_map(homicide, "2009-01-01", "2023-12-31", "murder_map_full_period.png")


# =====================================================================
# 7. WestSide_HotSpot.png / SouthSide_HotSpot.png  — 3×3 heatmap
# =====================================================================
def plot_heatmaps():
    print("\n7. Heatmap grids (WestSide / SouthSide)")
    chicago = gpd.read_file(os.path.join(DATA, "CommAreas_20250408.geojson"))
    chicago.loc[:, "community"] = chicago["community"].str.title()

    hom = pd.read_csv(os.path.join(CLEAN, "Homicide_six_hoods.csv"))
    ssp = pd.read_csv(os.path.join(CLEAN, "Shotspotter_six_hoods.csv"))
    st  = pd.read_csv(os.path.join(DATA, "transportation_20250415.csv"))

    hom['geometry'] = gpd.points_from_xy(hom['Longitude'], hom['Latitude'])
    ssp['geometry'] = gpd.points_from_xy(ssp['Longitude'], ssp['Latitude'])
    hom_gdf = gpd.GeoDataFrame(hom, geometry='geometry', crs='EPSG:4326')
    ssp_gdf = gpd.GeoDataFrame(ssp, geometry='geometry', crs='EPSG:4326')
    st['geometry'] = st['the_geom'].apply(wkt.loads)
    st_gdf = gpd.GeoDataFrame(st, geometry='geometry', crs='EPSG:4326')

    chicago = chicago.to_crs('EPSG:26971')
    hom_gdf = hom_gdf.to_crs('EPSG:26971')
    ssp_gdf = ssp_gdf.to_crs('EPSG:26971')
    st_gdf  = st_gdf.to_crs('EPSG:26971')

    def create_grid(geom, sz=150):
        x0, y0, x1, y1 = geom.bounds
        grids = []
        for x in np.arange(x0, x1, sz):
            for y in np.arange(y0, y1, sz):
                grids.append(Polygon([(x,y),(x+sz,y),(x+sz,y+sz),(x,y+sz)]))
        return gpd.GeoDataFrame({'geometry': grids}, crs='EPSG:26971')

    def assign_to_grid(pts, grid):
        grid = grid.copy()
        joined = gpd.sjoin(pts, grid, how='left', predicate='intersects')
        cnts = joined.groupby('index_right').size()
        grid['count'] = grid.index.map(cnts).fillna(0)
        return grid

    pre  = hom_gdf[hom_gdf['Date'] < '2017-01-13']
    post = hom_gdf[hom_gdf['Date'] >= '2017-01-13']

    def make_heatmap(hoods, filename):
        buf = 400

        # Pre-compute global maxes
        gmax_pre = gmax_post = gmax_ss = 0
        for h in hoods:
            ng = chicago[chicago['community'] == h].geometry.values[0]
            g1 = assign_to_grid(pre,     create_grid(ng))
            g2 = assign_to_grid(post,    create_grid(ng))
            g3 = assign_to_grid(ssp_gdf, create_grid(ng))
            gmax_pre  = max(gmax_pre,  g1['count'].max())
            gmax_post = max(gmax_post, g2['count'].max())
            gmax_ss   = max(gmax_ss,   g3['count'].max())

        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(16, 16))
        plt.subplots_adjust(wspace=0.08, hspace=0.18,
                            left=0.03, right=0.88, top=0.93, bottom=0.04)

        for idx, hood in enumerate(hoods):
            ax1, ax2, ax3 = axes[idx]
            ng = chicago[chicago['community'] == hood].geometry.values[0]
            nb = ng.buffer(buf)
            streets_clip = gpd.clip(st_gdf, gpd.GeoDataFrame(geometry=[ng], crs='EPSG:26971'))

            pg = assign_to_grid(pre,     create_grid(ng))
            og = assign_to_grid(post,    create_grid(ng))
            sg = assign_to_grid(ssp_gdf, create_grid(ng))

            def bg(ax):
                # contextily downloads basemap tiles over the network; if offline
                # the map still renders (without the street/terrain backdrop).
                try:
                    ctx.add_basemap(ax, crs='EPSG:26971',
                                    source=ctx.providers.CartoDB.Positron,
                                    alpha=0.45)
                except Exception as e:
                    print(f"    [warn] basemap tiles unavailable (offline?); "
                          f"rendering without basemap: {e}")
                ax.set_xlim(nb.bounds[0], nb.bounds[2])
                ax.set_ylim(nb.bounds[1], nb.bounds[3])
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)

            bg(ax1)
            pg.plot(ax=ax1, column='count', cmap='OrRd', legend=False,
                    alpha=0.9, vmin=0, vmax=gmax_pre)
            streets_clip.plot(ax=ax1, color='#444444', linewidth=0.3, alpha=0.35)
            gpd.GeoDataFrame(geometry=[ng], crs='EPSG:26971').boundary.plot(
                ax=ax1, color='black', linewidth=1.5)
            ax1.set_title(f'{hood}\nGun Homicide Pre-2017',
                          fontsize=11, fontweight='bold')

            bg(ax2)
            sg.plot(ax=ax2, column='count', cmap='Blues', legend=False,
                    alpha=0.9, vmin=0, vmax=gmax_ss)
            streets_clip.plot(ax=ax2, color='#444444', linewidth=0.3, alpha=0.35)
            gpd.GeoDataFrame(geometry=[ng], crs='EPSG:26971').boundary.plot(
                ax=ax2, color='black', linewidth=1.5)
            ax2.set_title(f'{hood}\nShotSpotter Alerts',
                          fontsize=11, fontweight='bold')

            bg(ax3)
            og.plot(ax=ax3, column='count', cmap='OrRd', legend=False,
                    alpha=0.9, vmin=0, vmax=gmax_post)
            streets_clip.plot(ax=ax3, color='#444444', linewidth=0.3, alpha=0.35)
            gpd.GeoDataFrame(geometry=[ng], crs='EPSG:26971').boundary.plot(
                ax=ax3, color='black', linewidth=1.5)
            ax3.set_title(f'{hood}\nGun Homicide Post-2017',
                          fontsize=11, fontweight='bold')

        # Colorbars
        sm_pre  = plt.cm.ScalarMappable(cmap='OrRd',
                    norm=plt.Normalize(0, gmax_pre))
        sm_ss   = plt.cm.ScalarMappable(cmap='Blues',
                    norm=plt.Normalize(0, gmax_ss))
        sm_post = plt.cm.ScalarMappable(cmap='OrRd',
                    norm=plt.Normalize(0, gmax_post))

        cb1 = fig.colorbar(sm_pre,  ax=axes[:, 0].tolist(), orientation='vertical',
                           fraction=0.03, pad=0.02)
        cb1.set_label('Gun Homicides Pre-2017', fontsize=10, weight='bold')
        cb1.ax.tick_params(labelsize=8)

        cb2 = fig.colorbar(sm_ss,   ax=axes[:, 1].tolist(), orientation='vertical',
                           fraction=0.03, pad=0.02)
        cb2.set_label('ShotSpotter Alerts', fontsize=10, weight='bold')
        cb2.ax.tick_params(labelsize=8)

        cb3 = fig.colorbar(sm_post, ax=axes[:, 2].tolist(), orientation='vertical',
                           fraction=0.03, pad=0.02)
        cb3.set_label('Gun Homicides Post-2017', fontsize=10, weight='bold')
        cb3.ax.tick_params(labelsize=8)

        side = "West Side" if "Austin" in hoods else "South Side"
        fig.suptitle(
            f'Gun Homicide & ShotSpotter Spatial Analysis — {side} Neighborhoods',
            fontsize=14, fontweight='bold', y=0.97)
        savefig(fig, filename)

    make_heatmap(['Austin', 'Humboldt Park', 'North Lawndale'],
                 'hotspot_west_side.png')
    make_heatmap(['Englewood', 'West Englewood', 'South Shore'],
                 'hotspot_south_side.png')


# =====================================================================
# MAIN
# =====================================================================
if __name__ == '__main__':
    plot_bar_chart()
    plot_monthly_timeseries()
    plot_monthly_mean()
    plot_correlation()
    plot_hourly()
    plot_murder_maps()
    plot_heatmaps()
    print("\nAll plots regenerated successfully.")
