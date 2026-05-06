#ifndef TIMING_HELPER_H
#define TIMING_HELPER_H

#include <iostream>
#include <chrono>
#include <sstream>
#include <vector>

namespace helper
{
    template <typename Clock = std::chrono::steady_clock>
    struct ChronoHelper
    {
        std::chrono::time_point<Clock> start;

        static std::chrono::time_point<Clock> now()
        {
            return Clock::now();
        }

        ChronoHelper()
        {
            reset();
        }

        void reset()
        {
            start = now();
        }

        template <typename Duration = std::chrono::milliseconds>
        Duration get_duration() const
        {
            return std::chrono::duration_cast<Duration>(now() - start);
        }

        std::size_t get_duration_milis() const
        {
            return get_duration<std::chrono::milliseconds>().count();
        }
    };

    template <typename Clock = std::chrono::steady_clock, typename Duration = std::chrono::milliseconds>
    struct ChronoProcess
    {
        ChronoHelper<Clock> total_chrono;
        ChronoHelper<Clock> next_chrono;
        Duration loop_duration;
        std::string dur_type = "ms";
        std::string speed_type = "speed";

        std::size_t nb_checks = 0;

        ChronoProcess(Duration loop_duration) : loop_duration(loop_duration)
        {
            reset();
        }

        ChronoProcess(std::size_t loop_duration_int) : loop_duration(loop_duration_int)
        {
            reset();
        }

        void reset()
        {
            total_chrono.reset();
            next_chrono.reset();
        }

        bool time_is_out()
        {
            if (next_chrono.template get_duration<Duration>() >= loop_duration)
            {
                next_chrono.reset();
                return true;
            }
            else
            {
                return false;
            }
        }

        // state is between 0 and 1
        void check_and_print(const float &state, const bool &erase_line = false)
        {
            nb_checks++;
            if (time_is_out() || state == 1.0)
            {
                if (erase_line)
                    std::cout << "\r";
                std::cout << "process : " << (int)(state * 100) << "% (" << total_chrono.template get_duration<Duration>().count() << " " << dur_type
                          << ", " << nb_checks << " " << speed_type << ")";
                if (erase_line)
                    std::cout << std::flush;
                else
                    std::cout << std::endl;
                nb_checks = 0;
            }
        }

        void print_end()
        {
            check_and_print(1.0, true);
            std::cout << std::endl;
        }
    };

    // defaults to blue
    // see https://en.wikipedia.org/wiki/ANSI_escape_code
    std::string cout_colors(std::vector<int> colors = {34})
    {
        if (colors.size() == 0)
            return "";
        std::stringstream s;
        s << "\033[" << colors[0];
        for (int i = 1; i < colors.size(); i++)
        {
            s << ";" << colors[i];
        }
        s << "m";
        return s.str();
    }
    std::string cout_color_reset()
    {
        return cout_colors({0});
    }

    namespace internal
    {

        struct Multi_ostream
        {
            std::vector<std::ostream *> ostreams;
            std::vector<bool> color_support;
            void add_trace_ostream(std::ostream &ostream, const bool does_support_color = true)
            {
                ostreams.push_back(&ostream);
                color_support.push_back(does_support_color);
            }
            void clear()
            {
                ostreams.clear();
                color_support.clear();
            }
            void print(const std::string &s) const
            {
                for (std::ostream *os : ostreams)
                {
                    (*os) << s;
                }
            }
            void print_command(std::ostream &(*command)(std::ostream &)) const
            {
                for (std::ostream *os : ostreams)
                {
                    (*os) << command;
                }
            }
            template <typename Command>
            void print_command(Command &command) const
            {
                for (std::ostream *os : ostreams)
                {
                    (*os) << command;
                }
            }
            void print_colors(std::vector<int> colors = {34}) const
            {
                for (std::size_t os_id = 0; os_id < ostreams.size(); os_id++)
                {
                    if (color_support[os_id])
                        (*(ostreams[os_id])) << cout_colors(colors);
                }
            }
            void print_color_reset() const
            {
                print_colors({0});
            }
        };

        struct Padding_info
        {
            Multi_ostream *ostream;
            std::size_t depth;
            std::string message;
            ChronoHelper<> chrono;

            Padding_info() {}
            Padding_info(Multi_ostream &ostream) : ostream(&ostream) {}
            Padding_info(Multi_ostream &ostream,
                         const std::size_t &depth, const std::string &message, const ChronoHelper<> &chrono)
                : ostream(&ostream), depth(depth), message(message), chrono(chrono)
            {
            }

            std::stringstream get_padding_stringstream() const
            {
                std::stringstream s;
                for (std::size_t i = 0; i < depth; i++)
                {
                    s << "    ";
                }
                return s;
            }

            void print_block_padding() const
            {
                ostream->print(get_padding_stringstream().str());
                ostream->print_command(std::flush);
            }
            void print_start() const
            {
                ostream->print(get_padding_stringstream().str());
                ostream->print_colors();
                std::stringstream s;
                s << "[" << message << "]";
                ostream->print(s.str());
                ostream->print_color_reset();
                ostream->print_command(std::flush);
            }

            void print_end_line(const std::string &end_message = "", const bool set_colors = true) const
            {
                if (set_colors)
                    ostream->print_colors();
                std::stringstream s;
                s << end_message << "(" << chrono.get_duration_milis() << " ms)";
                ostream->print(s.str());
                ostream->print_color_reset();
                ostream->print_command(std::endl);
            }
            void print_end_block() const
            {
                ostream->print(get_padding_stringstream().str());
                ostream->print_colors();
                std::stringstream s;
                s << "[" << message << "]";
                ostream->print(s.str());
                print_end_line("", false);
            }
        };

        struct Trace
        {
            std::vector<Padding_info> padding_infos;
            Multi_ostream ostreams;

            Trace()
            {
                ostreams.add_trace_ostream(std::cout);
            }

            std::string current_padding_string() const
            {
                return padding_infos.back().get_padding_stringstream().str();
            }

            void bLine(const std::string &message = "")
            {
                std::vector<Padding_info>::iterator last_padding =
                    padding_infos.emplace(padding_infos.end(),
                                          ostreams,
                                          padding_infos.size(),
                                          message,
                                          ChronoHelper<>());
                last_padding->chrono.reset();
                last_padding->print_start();
            }

            void bBlock(const std::string &message)
            {
                bLine(message);
                ostreams.print_command(std::endl);
            }

            void eLine(const std::string &message = "")
            {
                const internal::Padding_info &last_padding = padding_infos.back();
                last_padding.print_end_line(message);
                padding_infos.resize(last_padding.depth);
            }

            void eBlock()
            {
                const internal::Padding_info &last_padding = padding_infos.back();
                last_padding.print_end_block();
                padding_infos.resize(last_padding.depth);
            }

            void message(const std::string &message, const bool &new_line = true) const
            {
                ostreams.print(current_padding_string());
                ostreams.print_colors();
                ostreams.print(message);
                ostreams.print_color_reset();
                ostreams.print_command(std::flush);
                if (new_line)
                    ostreams.print_command(std::endl);
                else
                    ostreams.print_command(std::flush);
            }
        };
        static Trace global_trace;

    } // namespace internal

    std::string current_padding_string()
    {
#ifdef VERBOSE
        return internal::global_trace.current_padding_string();
#endif
        return "";
    }

    void bLine(const std::string &message = "")
    {
#ifdef VERBOSE
        internal::global_trace.bLine(message);
#endif
    }

    void eLine(const std::string &message = "")
    {
#ifdef VERBOSE
        internal::global_trace.eLine(message);
#endif
    }

    void bBlock(const std::string &message)
    {
#ifdef VERBOSE
        internal::global_trace.bBlock(message);
#endif
    }

    void eBlock()
    {
#ifdef VERBOSE
        internal::global_trace.eBlock();
#endif
    }

    void message(const std::string &message, const bool &new_line = true)
    {
#ifdef VERBOSE
        internal::global_trace.message(message, new_line);
#endif
    }

    void clear_trace_ostream()
    {
#ifdef VERBOSE
        internal::global_trace.ostreams.clear();
#endif
    }

    void add_trace_ostream(std::ostream &ostream, const bool does_support_color = true)
    {
#ifdef VERBOSE
        internal::global_trace.ostreams.add_trace_ostream(ostream, does_support_color);
#endif
    }
} // namespace helper

#endif // TIMING_HELPER_H