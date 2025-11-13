        # if plotGXG:
        #     ax = axs[0]
        #     hyears = self.GXG['hyears']
            
        #     ax.set_title(ax.get_title() + f" {hyears[0]}-{hyears[-1]}")
            
        #     for ip in parcels:
        #         ax.plot(self.tdata, self.avgHds[:, 0, ip], label=f"parcel[{ip}]")

        #         ax.axhline(self.GXG['GHG'][ip], c='b', label='GHG')
        #         ax.axhline(self.GXG['GVG'][ip], c='g', label='GVG')
        #         ax.axhline(self.GXG['GLG'][ip], c='r', label='GLG')
                
        #         for iy, hyear in enumerate(self.VG3):
        #             if iy == 0:
        #                 labels = ["HG3", "VG3", "LG3"]
        #             else:
        #                 labels = ["", "", ""]
        #             ax.plot(self.HG3[hyear]['t'][:, ip], self.HG3[hyear]['h'][:, ip], 'b^', label=labels[0])                    
        #             ax.plot(self.VG3[hyear]['t'][:, ip], self.VG3[hyear]['h'][:, ip], 'go', label=labels[1])                    
        #             ax.plot(self.LG3[hyear]['t'][:, ip], self.LG3[hyear]['h'][:, ip], 'rv', label=labels[2])
                    
        #     ax.legend()
            
        #     plot_hydrological_year_boundaries(ax=ax, tindex=self.tdata.index)
